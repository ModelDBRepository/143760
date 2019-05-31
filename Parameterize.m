function [Alpha, AlphaApprox, TauKappa, Beta, Kappa, TauJ] = Parameterize(RelativeSpread, Chronaxie, TauSum, Threshold, Jitter, KnownInput)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameterize.m
% Joshua Goldwyn
% 5/10/11
%  
% This function parameterizes the AN model presented in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Values used in simulations
%[Alpha, AlphaApprox, TauKappa, Beta, Kappa, TauJ] = Parameterize(.0487, 276, 250, .852, 85.5,[0 0 0 0 0]); 
%[Alpha, AlphaApprox, TauKappa, Beta, TauJ] =   [25.6337   24.5196 325.3700 0.3330   9.3649 96.9336]
%%%


% Some needed variables
dt = 1; % mu sec, for evaluating integrals and filters
MaxPhaseDur = 2000;  % Used in Chronaxie mu sec

% Pulses for Summation Experiment (Cartee et al 2006)
SummationPhaseDur = 50;
SummationPulseShape = 'pseudo';
IPIall = [100 200 300];

% Pulses for Threshold Experiment (Miller et al 2001 biphasic experiments)
ThresholdPhaseDur = 40;
ThresholdPulseShape = 'bi';

% Pulses for Jitter Experiment (Miller et al 2001 biphasic experiments)
JitterPhaseDur = 40;
JitterPulseShape = 'bi';

% Approximation of Alpha vs RS
AlphaApproximation = @(rs) rs.^(-1.0587);
% %%%% NOTE THIS CODE GIVES THE APPROXIMATION:
% RSfunc = @(a) sqrt(gamma(1+2./a)./(gamma(1+1./a)).^2 - 1);
% Alpha = linspace(1,100,1E3);
% p = fminsearch(@(p) norm(RSfunc(Alpha).^p-Alpha),-1);
% %%%%%%%%%%%%%%%%


if KnownInput(1)>0
    Alpha = KnownInput(1);
else
%%% First Compute Alpha for given RS
% RS as a function of alpha using, using st dev / mean of Weibull distribution
RSfunc = @(a) sqrt(gamma(1+2./a)./(gamma(1+1./a)).^2 - 1);
Alpha = fminsearch(@(a) norm(RSfunc(a) - RelativeSpread), 10);
end
AlphaApprox = AlphaApproximation(RelativeSpread);


if KnownInput(2)>0
    TauKappa = KnownInput(2);
else
%%% Next Compute TauKappa for given Chronaxie
% Maximum Phase Duration Used in Van den Honert and Stypulkowski 1984
TauKappa = fminsearch(@(tk) ChronaxieFit(tk, Chronaxie, MaxPhaseDur, AlphaApprox, dt), 100);
end

if KnownInput(3)>0
    Beta = KnownInput(3);
else
%%% Next Compute Beta for given Summation Time Constant
Beta = fminbnd(@(b) SummationFit(b, TauSum, SummationPhaseDur, AlphaApprox, TauKappa, dt, SummationPulseShape, IPIall),0,1);
end

if KnownInput(4)>0
    Kappa = KnownInput(4);
else
%%% Next Compute Kappa for given Threshold
tend = 5000;
t = 0:dt:tend;
D = ThresholdPhaseDur;  %phase duration
switch ThresholdPulseShape
    case 'bi'
        w = zeros(1,length(t)); w(1:D/dt) = 1; w(D/dt+1:2*D/dt)=-Beta; %biphasic
end

wFold = 0;
wold = 0;
for i=1:length(t)
    wFilter(i) = wFold*exp(-dt/TauKappa)+(dt / (2*TauKappa))*(wold*exp(-dt/TauKappa)+w(i));%conv(1/TauKappa * exp(-t/TauKappa), w)*dt;
    wFold = wFilter(i);
    wold = w(i);
end
W = wFilter(1:length(t)).*(wFilter(1:length(t))>0); % output of stimulus filter, no kappa
intW = trapz(W.^AlphaApprox)*dt; % integrate response

Kappa = (log(2)/intW).^(1/AlphaApprox)./ Threshold;


end


if KnownInput(5)>0
    TauJ = KnownInput(5);
else
%%% Last, compute TauJ for Jitter Filter
TauJ = fminsearch(@(tj) JitterFit(Jitter, tj, JitterPhaseDur, TauKappa, Beta, AlphaApprox, Kappa, Threshold,dt), 100);
end

end



function [CAerror] = ChronaxieFit(TauKappa, Chronaxie, MaxPhaseDur, Alpha, dt)
    % TauKappa will be fit based on this function

    % Total time to evaluate FE
    tend = 5000; 
    
    % Phase durations
    D0 = Chronaxie/dt;
    Dmax = MaxPhaseDur/dt;
    
    t = 0:dt:tend;

    % waveform functions, always monophasic
    w1 = zeros(1,length(t)); w1(1:D0/dt) = 1; 
    w2 = zeros(1,length(t)); w2(1:Dmax/dt) = 1; 
  
    % Apply K filter
    w1Filter = conv(1/TauKappa * exp(-t/TauKappa), w1)*dt;
    w2Filter = conv(1/TauKappa * exp(-t/TauKappa), w2)*dt;
    W1 = w1Filter(1:length(t)).*(w1Filter(1:length(t))>0); % output of stimulus filter, no kappa because cancels out later
    W2 = w2Filter(1:length(t)).*(w2Filter(1:length(t))>0); % output of stimulus filter, no kappa because cancels out later

    % Apply nonlinearity and Integrate Responses
    intW1 = trapz(W1.^Alpha)*dt;
    intW2 = trapz(W2.^Alpha)*dt;
    
    % Error function that can be used to determine TauKappa
    CAerror = norm(2^Alpha - (intW2/intW1));
    
end

function [SumError] = SummationFit(Beta, TauSum, PhaseDur, Alpha, TauKappa, dt, StimulusType, IPIall)

    % Total time to evaluate FE
    tend = 5000; 

    % Pulse Tain Parameters
    D =  PhaseDur;

    % t
    t = 0:dt:tend;
 
for i=1:length(IPIall)
    IPI = IPIall(i);
    
    % PseudoMonophasic
    switch StimulusType
        case 'pseudo'
            tau = 5E6; % Cartee time constant defined by hardware
            w1 = zeros(1,length(t)); w1(1:D/dt) = 1; 
            tt = t(D/dt+1:IPI/dt) - D;
            w1(D/dt+1:IPI/dt)= w1(D/dt+1:IPI/dt) -Beta * D * exp(-(tt-D)/tau) / (tau*(1-exp(-(IPI-D)/tau)));
            w2= w1; w2(IPI/dt+1:IPI/dt+D/dt) = w2(IPI/dt+1:IPI/dt+D/dt)+1;  % 2nd pulse
            w2((IPI+D)/dt+1:(2*IPI)/dt)= w2((IPI+D)/dt+1:(2*IPI)/dt) - Beta*D*exp(-(tt-D)/tau) / (tau*(1-exp(-(IPI-D)/tau)));

        case 'bi' % Biphasic
            w1 = zeros(1,length(t)); w1(1:D/dt) = 1; 
            w1(D/dt+1:(2*D)/dt)= w1(D/dt+1:(2*D)/dt) - Beta;
            w2= w1; w2(IPI/dt+1:IPI/dt+D/dt) = w2(IPI/dt+1:IPI/dt+D/dt)+1;  % 2nd pulse
            w2((IPI+D)/dt+1:(IPI+2*D)/dt) = w2((IPI+D)/dt+1:(IPI+2*D)/dt)-Beta;  % 2nd pulse
    end

        % Apply K filter
        w1Filter = conv(1/TauKappa * exp(-t/TauKappa), w1)*dt;
        w2Filter = conv(1/TauKappa * exp(-t/TauKappa), w2)*dt;
        W1 = w1Filter(1:length(t)).*(w1Filter(1:length(t))>0); % output of stimulus filter, no kappa
        W2 = w2Filter(1:length(t)).*(w2Filter(1:length(t))>0); % output of stimulus filter, no kapp

        intW1 = trapz(W1.^Alpha)*dt;
        intW2 = trapz(W2.^Alpha)*dt;

        ThreshRatio(i) = (intW1/intW2)^(1/Alpha);
    end

    % Error Using Equation in Cartee et al
    SumError = norm( (1-.5*exp(-IPIall/TauSum)) - ThreshRatio);

end





function [JitterError] = JitterFit(jit, tauj, PhaseDur, TauKappaFit, BetaFit, AlphaFit, KappaFit, ThreshVal,dt)

    tend = 5000;
    t = 0:dt:tend;
    D = PhaseDur;  %phase duration
    w = zeros(1,length(t)); w(1:D/dt) = 1; w(D/dt+1:2*D/dt)=-BetaFit; %biphasic
    wFilter = conv(1/TauKappaFit * exp(-t/TauKappaFit), w)*dt;
    W = wFilter(1:length(t)).*(wFilter(1:length(t))>0); % output of stimulus filter, no kappa
 
    JWFilter = conv(1/tauj* exp(-t/tauj), W.^AlphaFit)*dt;
    JW = JWFilter(1:length(t)).*(JWFilter(1:length(t))>0);
    lambda = (KappaFit*ThreshVal)^AlphaFit * JW;
    Lambda = cumtrapz(lambda)*dt;
    f = 2*lambda.*exp(-Lambda);
    jj = sqrt(trapz(t.^2.*f) - trapz(t.*f).^2 );

    JitterError = norm(jj -jit);

end