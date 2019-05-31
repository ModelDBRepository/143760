function I = Current(t,P)

% Define the current input for model

% t is time
% P stores any needed parameter values

PulseLevel = P(1); 
PhaseDuration = P(2); 
PulseRate = P(3);

%%% Here is an example for biphasic pulse trains
I = zeros(size(t));
I( mod(t,1E6/PulseRate)<2*PhaseDuration ) = -PulseLevel;
I( mod(t,1E6/PulseRate)<PhaseDuration ) =    PulseLevel;
