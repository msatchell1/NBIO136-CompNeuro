% Tutorial_8_2.m
%
%  Models synaptic cometition between groups of inputs,
%  involved in the formation of ocular dominance stripes.
%  The inputs are grouped into two, with each group having
%  a time-varying (sinusoidal) firing rate. The rates of the two groups are
%  phase-offset. Each input is an inhomogeneous Poisson spike train whose
%  rate is that of its corresponding group. Each spike produces a step
%  change and exponential decay of an excitatory conductance in a 
%  leaky-integrate-and-fire cell, which fires spikes.
%
%  Synaptic strengths between inputs and the postsynaptic LIF neuron
%  are updated according to a rule for spike-timing dependent plasticity.
%
%  This code produces the results for Tutorial 8.2 in the book:
%  An Introductory Course in Computational Neuroscience,
%  by Paul Miller.
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

question_flag = 5;      % Questions 5 and 6 require altered inputs.
Ntrials = 1;          % number of runs with plasticity between each run

tmax = 0.5;             % time for a run (sec)
dt = 0.0001;            % timestep for integration (sec)
t = 0:dt:tmax;          % time vector (time in sec)

%% Set up the inputs to the neuron 
Ninputs = 50;           % total number of input cells
phase_offset = zeros(1,Ninputs);  % phase of cycle for each input, default 0
freq = 20;              % Default frequency of oscillating input (Hz)

switch question_flag
    case 5
        phase_offset(Ninputs/2+1:Ninputs) = pi;
    case 6
        freq = 4;              % reduced frequency of sinusoid input (Hz)
        N_correlated =Ninputs/5;
        % Random phase offset for the uncorrelated inputs
        phase_offset(N_correlated+1:Ninputs)  = 2*pi*rand(1,Ninputs-N_correlated);
    otherwise
        phase_offset(Ninputs/2+1:Ninputs) = pi;
end

tausyn = 0.002;         % synaptic timeconstant (sec)
itau = ceil(tausyn/dt); % number of timesteps in tau

rapp = zeros(length(t),2);  % input spike rate as function of time for each of two groups
rmax = 60;                  % peak input spike rate

%% Set up initial synaptic input strengths
W0 = 0.5e-9;                 % initial average synaptic strengths
Wrand = 0.025e-9;
% Add 5% random variation to the average
W = W0*ones(Ninputs,1) + Wrand*randn(Ninputs,1); 

%% Set up rules for changes in synaptic strengths
dW_stdpp = 0.02e-9;     % maximum change in strength for each pairing (S)
dW_stdpm = 0.025e-9;    % maximum change in strength for each pairing (S)
tau_stdpp = 0.020;      % time window of 2 spikes for potentiation (sec)
tau_stdpm = 0.020;      % time window of 2 spikes for depression (sec)
Wmax = 4*W0;            % maximum allowed synaptic strength
W = min(W,Wmax);        % ensure W is no higher than Wmax
W = max(W,0);           % and no lower than zero

%% Parameters for the LIF neuron being simulated

Vth = -0.050;           % threshold voltage (V)
E_L = -0.070;           % leak voltage (V)
E_syn = 0.000;          % Reversal potential of synaptic inputs (V)
Vreset = -0.080;        % reset voltage (V)
Cm = 100e-12;           % specific membrane capacitance (F)
G_L = 5e-9;             % specific leak conductance (S)

%% Begin the loop through all trials with plasticity implemented by the 
%  block method at the end of each trial

for trial = 1:Ntrials   % begin trials
    
    %% Produce the two sets of inputs with an offset that alternates across trials.
    %  Inputs are excitatory conductances due to Poisson inputs with
    %  sinusoidally vary input rates. The two sinusoids correspond to the
    %  two eyes with phase lag that can be altered. Default is out of phase
    %  by pi, in which case "leading" or a "lagging" are the same.
    
    for i = 1:Ninputs
        if mod(trial,2) == 0                        % on even trials
            rapp(:,i) = 0.5*rmax*(1+sin(2*pi*freq*t + phase_offset(i)));   % input rate for group 1
        else;                                       % on odd trials
            rapp(:,i) = 0.5*rmax*(1+sin(2*pi*freq*t - phase_offset(i))); % follow those of group 1            
        end
    end
    
    spikesin = zeros(Ninputs,length(t));    % series of spikes for each input
    spikesout = zeros(size(t));             % series of spikes generated by LIF neuron
    G_syn = zeros(size(t));                  % input current to LIF neuron
    
    for i = 1:Ninputs                     % first group of inputs
        for j = 1:length(t)
            if rapp(j,i)*dt > rand(1)       % random Poisson process at rate Iapp
                spikesin(i,j) = 1;          % to generate spikes at inputs
                
                kmax = min(j+5*itau,length(t));     % how long to add current
                for k = j+1:kmax                    % from the spike
                    G_syn(k) = G_syn(k) + ...         % add current from input i
                        W(i)*exp(-dt*(k-j)/tausyn);
                end
            end
        end
    end
    Gmean = mean(G_syn)
    % Completed generation of set of input current from input spike trains
    %% Now see how the LIF neuron responds to these inputs
    
    V = zeros(size(t));     % membrane potential of LIF neuron
    V(1) = E_L;              % initialize at leakage potential
    
    for i = 2:length(t)     % integrate through time
        
        Gtot = G_L + G_syn(i-1);
        Vinf = (E_L*G_L + E_syn*G_syn(i-1))/Gtot;           % calculate steady state voltage
        V(i) = Vinf + (V(i-1)-Vinf)*exp(-dt*Gtot/Cm);  % and integrate towards it
        
        if V(i) >= Vth                  % if voltage passes threshold
            spikesout(i) = 1;           % generate spike in LIF neuron
            V(i) = Vreset;              % and reset the membrane potential
        end
    end
    
    s1 = find(spikesout);   % series of spike times of LIF cell
    N1 = length(s1)         % total number of LIF spikes
    for cell2 = 1:Ninputs;  % for each input
        dWsum = 0.0;        % initialize change in weight as zero
        
        s2 = find(spikesin(cell2,:));   % series of spike times for input cell
        N2 = length(s2);                % number of input spikes from that cell
        
        for i = 1:N1                    % i runs through LIF spikes
            for j = 1:N2                % j runs through input cell spikes
                deltat = dt*(s1(i)-s2(j)); % difference in spike times
                if (deltat > 0 )        % if pre before post
                    dWsum = dWsum + dW_stdpp*exp(-deltat/tau_stdpp); % do LTP with decay window
                end
                if ( deltat < 0 )       % if post before pre
                    dWsum = dWsum - dW_stdpm*exp(deltat/tau_stdpm); % reduce synaptic weight
                end
            end
        end
        
        W(cell2) = W(cell2) + dWsum;    % now change the synaptic strength of that cell
    end
    
    W = min(W,Wmax);                        % now ensure no W is above maximum allowed
    W = max(W,0);                           % and no W is less than zero
    
    switch question_flag
        case 6
            input1 = mean(W(1:N_correlated));
            input2 = mean(W(N_correlated+1:end));
        otherwise
            input1 = mean(W(1:Ninputs/2));           % input1 is average of weights from group 1
            input2 = mean(W(Ninputs/2+1:Ninputs));   % input 2 is average of weights form group 2
    end
    
    strengths(trial,:) = [input1 input2];
    
    if ( trial == 1 ) || ( mod(trial,50) == 0 )
        figure(1)
        subplot(2,1,1)
        plot(t,V);              % plot membrane potential as function of time
        ylabel('Membrane potential (V)')
        subplot(2,1,2)
        plot(t,G_syn)
        xlabel('Time (sec)')
        ylabel('G_{syn} (Total)')
        drawnow
        
        figure(2)
        plot(W,'o')                         % plot the set of wights in Fig.2
        axis([0 Ninputs 0 Wmax])
        xlabel('Input number')
        ylabel('Synaptic strength (S)')
        drawnow
        
        figure(3)
        plot(strengths)
        xlabel('Trial Number')
        ylabel('Mean strength (S)')
        axis([0 Ntrials 0 Wmax]);
        legend('Group 1', 'Group 2')
        drawnow
                
        figure(5)
        plot(phase_offset,W,'o')
        xlabel('Phase Offset')
        ylabel('Synaptic Strength')
        drawnow
 
    end
end