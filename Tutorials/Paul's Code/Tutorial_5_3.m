% Tutorial_5_3.m
% This code has 2 leaky integrate and fire neuron with
% added noise term scaled by sigma.
% Neurons are coupled by synapses with simple 
% exponential decay time constant.
% The circuit without depression is bistable, with either cell active while
% the other is inactive. Noise or depression cause transitions (irregular
% or regular) between the two states.
%
% This code is a solution of Tutorial 5.3 in the textbook 
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

noise_flag = 1;     % set to zero to remove noise
dep_flag = 1;       % set to zero to remove synaptic depression
distribution_flag = 1;  % set to 1 to produce a distribution of state durations


E=-70e-3;           % Leak potential
Rm=10e6;            % Membrane resistance
C_m = 1e-9;         % Membrane capacitance

Vth=-54e-3;         % Threshold
Vreset=-80e-3;      % Reset voltage
Iappmax = 3e-9;     % Maximum value of applied current

Esyn12 = -70e-3;         % Reversal potential for synapse
Esyn21 = -70e-3;         % Use 0 for excitatory and -70e-3 for inhibitory

g21 = 1e-6;         % Conductance of synapse from 2 to 1
g12 = 1e-6;         % Conductance of synapse from 1 to 2

tau_s1 = 0.010;     % time constant of synapse from 1 to 2
tau_s2 = 0.010;     % time constant of synapse from 2 to 1

if (dep_flag )
    pr = 0.2;           % release probability (used if dep_flag is 1)
else
    pr = 1;
end
tau_d = 0.2;       % depression time constant (used if dep_flag is 1)

if ( noise_flag)
    if ( dep_flag )
        sigma = 5e-12;
    else
        sigma = 5e-11;       % sigma scales the noise term in the current
    end
else
    sigma = 0;
end

dt=0.0001;          % Time step of 0.1ms is a little high
tmax=6;             % Maximum time used to calculate rate
tpulse=3.0;         % This starts the pulse at the beginning of the trial
lengthpulse=0.1;    % and finishes it at the end of the trial.

if (distribution_flag ) % If a distribution of state durations is required
    tmax = 2000;        % Simulate for much longer 
end

Nt=tmax/dt +1;      % Nt is number of time points
T=0:dt:tmax;        % Time vector

Iapp = 3e-9;        % Extra current used to switch states
if ( noise_flag || distribution_flag) % If noise is present OR a state-durations are required
    Iapp = 0;       % Do not use a current pulse to force transitions
end

%% Initialize variables as separate vectors, one for each neuron
Iapp0 = 2e-9; % Baseline current to both neurons
Iapp1=zeros(size(T));       % Applied current to neuron 1
Iapp2=zeros(size(T));       % Applied current to neuron 2
I1=zeros(size(T));          % Total current to neuron 1
I2=zeros(size(T));          % Total current to neuron 2
V1=zeros(size(T));          % Membrane potential of neuron 1
V2=zeros(size(T));          % Membrane potential of neuron 2
D1=ones(size(T));   % initialize depression variable as 1
D2=ones(size(T));   % initialize depression variable as 1
s1=zeros(size(T));          % Synaptic gating variable of neuron 1
s2=zeros(size(T));          % Synaptic gating variable of neuron 2

noise = randn(2,length(T))*sigma/sqrt(dt);  % 2-rows of noise to add to I

V1(1)=E;                    % Begin at rest (leak) voltage
V2(1)=E;                    % Begin at rest (leak) voltage
tref=0.002;                 % Refractory period
lastspike1=-2*tref;         % Begin outside the last refractory period
lastspike2=-2*tref;         % Begin outside the last refractory period
spikes1=zeros(size(T));     % Binary vector of spikes for neuron 1
spikes2=zeros(size(T));     % Binary vector of spikes for neuron 2

%% Set the applied currents to be nonzero during a pulse
for i = 1 : round(lengthpulse/dt)
    Iapp1(i) = Iapp;     % Set initial current pulse to neuron 1
end;

for i = round(tpulse/dt)+1 : round((tpulse+lengthpulse)/dt)
    Iapp2(i) = Iapp;     % Set second current pulse to neuron 2
end;

%% Now integrate through time
for i = 2:Nt;    
    
    I1(i-1) = Iapp0+Iapp1(i-1) + g21*s2(i-1)*(Esyn21-V1(i-1)) + noise(1,i-1);
    V1(i) = V1(i-1) + dt*( (E-V1(i-1))/Rm + I1(i-1) )/C_m;

    I2(i-1) = Iapp0+Iapp2(i-1) + g12*s1(i-1)*(Esyn12-V2(i-1)) + noise(2,i-1);
    V2(i) = V2(i-1) + dt*( (E-V2(i-1))/Rm + I2(i-1) )/C_m;

    D1(i) = D1(i-1) + (1-D1(i-1))*dt/tau_d;     % Forward Euler update of D1
    D2(i) = D2(i-1) + (1-D2(i-1))*dt/tau_d;     % Forward Euler update of D2
    
    s1(i) = s1(i-1)*exp(-dt/tau_s1);    % exponential Euler for decay of s1
    s2(i) = s2(i-1)*exp(-dt/tau_s2);    % exponential Euler for decay of s2
    
    if ( T(i) < lastspike1 + tref ) % is still in the refractory period
        V1(i) = Vreset;             % keep at reset
    end
    if V1(i) > Vth       % if voltage is above threshold
        V1(i) = Vreset;  % set to reset
        spikes1(i) = 1;  % record a spike
        lastspike1=T(i); % time of last spike
        s1(i) = s1(i)+pr*D1(i)*(1-s1(i));   % Increase synaptic gating variable
        if ( dep_flag)
            D1(i) = D1(i)*(1-pr);           % Reduce depression variable
        end
    end
    
    if ( T(i) < lastspike2 + tref ) % is still in the refractory period
        V2(i) = Vreset;             % keep at reset
    end
    if V2(i) > Vth       % if voltage is above threshold
        V2(i) = Vreset;  % set to reset
        spikes2(i) = 1;  % record a spike
        lastspike2=T(i); % time of last spike
        s2(i) = s2(i)+pr*D2(i)*(1-s2(i));   % Increase synaptic gating variable
        if ( dep_flag)
            D2(i) = D2(i)*(1-pr);           % Reduce depression variable
        end
    end
end;    % end of time loop, for i = 2:Nt


%% Plot the results
% First set up font widths etc/
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

% Figure 1 has membrane potentials of the two cells superposed
figure(1)
clf
plot(T,V1)
hold on
plot(T,V2,'r')
axis([0 4 -0.08 -0.02])
% Figure 2 has two panels, with spike times of each cell in top panel and
% the synaptic gating variable of each cell in the bottom panel
figure(2)
clf
subplot(2,1,1)
plot(T,spikes1)
axis([0 4 0 1])
hold on
plot(T,spikes2,'r')
legend('spikes1', 'spikes2')
subplot(2,1,2)
axis([0 4 0 1])
hold on
plot(T,s1)
plot(T,s2,'r')
legend('s1', 's2')
xlabel('Time (sec)')

%% Calculate the state-transition times and distribution of durations in each state
start_testing = 0.5;    % Ignore initial transient of 0.5s
if (distribution_flag)  % If distribution of durations is required
    num_transitions = 0;    % Initialize counter for number of state transitions
    transition_times = [];  % Initialize record of transition times 
    
    states = zeros(size(T));    % Will record 1s and 2s for current state
    states(s1>s2) = 1;
    states(s2>s1) = 2;
    
    for i = round(start_testing/dt):Nt      % Loop through time-points
        if ( states(i) ~= states(i-1) )     % if state changes at time-point
            num_transitions = num_transitions + 1;      % add a transition
            transition_times(num_transitions) = T(i);   % record the time
        end
    end
    
    state_durations = diff(transition_times);   % Set of durations of states
    
    % In this example, the system is symmetric so transitions from 1 to 2
    % are equivalent to transitions from 2 to 1 and can be combined. 
    % In general, odd state durations and even state durations should be
    % separated.
    figure(3)
    histogram(state_durations,31)
end
