% Sample code to generate a spike train
% Point process model developed by Goldwyn, Rubinstein, Shea-Brown for
% response of auditory nerve fiber to cochlear implant stimulation
% See Goldwyn, Rubinstein, Shea-Brown "A point process framework for modeling electrical stimulation of the auditory nerve" arXiv:1201.5428
% Last updated: June 2012 (JHG)

% User can define:
% Neural Parameters: threshold, relative spread, chronaxie, summation time constant, jitter, refractory effects on threshold and relative spread
% Simulation Parameters: Length of stimulus, pulse rate, phase duration, current level per pulse

% Spike train Output variables:
%  SpikeCount (number of spikes) and SpikeTrain (list of spike times in  micro sec)

close all
clear all

% Neural parameters defining the model
RelativeSpread = 0.0487;
Chronaxie = 276; % micro sec
TauSum = 250; % micro sec
Threshold = 0.852; % mA
Jitter = 85.5;  % micro sec
AbsRef = 332; % micro sec
RelRef = 411; % Time scale of relative refractory period (micro sec)
AbsRS = 199. ;
RelRS = 423. ;
ThresholdPhaseDuration = 40;  % Phase duration used for stimulation that defines threshold

% Compute associated model parameters
[Alpha, AlphaApprox, TauKappa, Beta, Kappa, TauJ] = Parameterize(RelativeSpread, Chronaxie, TauSum, Threshold, Jitter,[0 0 0 0 0]);

% Define Refractory Functions (time in these functions is in micro s)
theta_func = @(theta,t) 99999.*(t<AbsRef) + (t>AbsRef)*theta / (1.0 - exp(-(t-AbsRef)/RelRef)) ;
rs_func = @(rs,t) min(.5, rs*(t<AbsRef) + (t>AbsRef)*rs / ( 1.0 - exp(-(t-AbsRS)/RelRS)));

% Build look up table that will be used to determine kappa for varying values of alpha
intW = zeros(5000,1);
tend = 5000;
dt = 1.;
tt = 0:dt:tend;
for i=1:size(intW)
    ialpha = i/100.;
    w = zeros(1,length(tt)); 
    w(1:ThresholdPhaseDuration/dt) = 1; 
    w(ThresholdPhaseDuration/dt+1:2*ThresholdPhaseDuration/dt)=-Beta; %biphasic
    wFold = 0;
    wold = 0;
    for ii=1:length(tt)
        wFilter(ii) = wFold*exp(-dt/TauKappa)+(dt / (2*TauKappa))*(wold*exp(-dt/TauKappa)+w(ii));
        wFold = wFilter(ii);
        wold = w(ii);
    end
    W = wFilter(1:length(tt)).*(wFilter(1:length(tt))>0); % output of stimulus filter, no kappa
    intW(i) = trapz(W.^ialpha)*dt; % integrate response
end

%%  This double-% symbol indicates that the script can be run in two pieces by pressing command+return with cursor in the desired box
% The first half other script is slow, but only needs to be run once per "neuron"
% The second half (everything below here) runs the spike train generator

tic;  % Start timer 

% Simulation time (micro sec)
t_begin = 0;
t_end = 1E6; 
dt = 1.;
t = t_begin:dt:t_end;
nt = length(t);

% Stimulus (Here is an example of a train of constant current level, biphasic pulses)
CurrentLevel = 0.462; % mA 
PhaseDuration = 40;   %micro  sec
PulseRate = 5000;  % pulses per second
P  = [CurrentLevel PhaseDuration PulseRate]; % Vector of Parameter Values passed in to Current.m
I = Current(t,P);  % Using Current.m to define stimulus


%%%%%%%%%%%%% Run Point Process Model And Record Spike Times %%%%%%%%%%%%%
% Initial Values
v = 0;
w = 0;
ci =0;
Integrate_ci = 0;
TimeSinceSpike = 1E6;  % Make it big if want no spike history at stimulus onset
AlphaVal = AlphaApprox;
ThresholdUpdate = Threshold;
RelativeSpreadUpdate= RelativeSpread;
SpikeCount = 0;
SpikeTrain = zeros(ceil(t_end/PulseRate),1);
r = -log(rand);
Iin = 0;

for i=2:nt
    
  TimeSinceSpike = TimeSinceSpike+dt;
  
  Iin_old = Iin;
  if (I(i-1)>0) % Positive phase
      Iin = I(i-1);
  elseif (I(i-1)<0) % Negative phase
      Iin = Beta*I(i-1);
  else
      Iin = 0;
  end
  
  if (TimeSinceSpike < AbsRef) % in absolute refractory period
    v = 0;
  else % Update state variable 
            v = v*exp(-dt/TauKappa) + (dt*Kappa/(2.*TauKappa)) * (Iin_old*exp(-dt/TauKappa) + Iin); % Stimulus filter, convolution using Trapezoid method
  end
  
  % Apply nonlinearity
  w_old = w;  % Save this old value for later calculation
  w = max(0,v)^AlphaVal;

  % Apply Jitter Filter
  ci_old = ci; % Conditional Intensity Value
  ci = ci*exp(-dt/TauJ) + dt /(2.*TauJ) * (w + w_old*exp(-dt/TauJ));
  
  % Integrate Conditional Intensity with Trapezoid method
  Integrate_ci = Integrate_ci + dt*ci;%(dt *(ci+ci_old))/2; 
  
  % Check for spikes
  if (Integrate_ci > r)  % SPIKE!

      % Record Spike
      SpikeCount = SpikeCount + 1;
      SpikeTrain(SpikeCount) = t(i);

      % Reset dynamical variables
      v=0; w=0; ci=0; Integrate_ci=0;  
      TimeSinceSpike = 0;
      
      % Draw New Random Number
      r = -log(rand);

  end
  
  % Update Spike Dependent Parameter Values at onset of next pulse
  if mod(t(i),1E6/(PulseRate*dt))==0 % new pulse

      RelativeSpreadUpdate = rs_func(RelativeSpread,TimeSinceSpike);
      ThresholdUpdate = theta_func(Threshold,TimeSinceSpike);
      AlphaVal = RelativeSpreadUpdate^(-1.0587); % approximation to exact expression       
      if TimeSinceSpike>AbsRef
          Kappa = (log(2)/intW(round(AlphaVal*100))).^(1/AlphaVal)./ ThresholdUpdate;
      else
          Kappa = 0;
      end
      
  end
  
end  % End loop over time steps


% Plot ISI histogram
maxt = 25;
[y,x] = hist(diff(.001*SpikeTrain(1:SpikeCount)),0:.05:maxt);
set(gca,'FontSize',18)
plot(x,y,'-','LineWidth',1)
xlabel('Interspike interval (ms)')
ylabel('Number of occurrences')
title('ISI histogram', 'FontSize',20)
xlim([0 maxt])

display(['Run Time = ', num2str(toc)])
display(['Spike/Second = ',num2str(SpikeCount/t_end*1E6)])