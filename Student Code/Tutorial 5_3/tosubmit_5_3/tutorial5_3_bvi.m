%tutorial 5.3 part a iv
%Sayaka (Saya) Minegishi
% minegishis@brandeis.edu
% Mar 12 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%% Simulation parameters (identical across cells)
dt = 0.0001;                % time-step
tmax = 13; %maximum simulation time (s)
t = 0:dt:tmax;               % vector of time-points (s)


istart1 = 0;          % time applied current starts (sec) for cell1
istart2 = tmax/2;          % time applied current starts (sec) for cell 2
ilength=0;    % length of applied current pulse (sec)

Ie1= 3e-9;          % magnitude of applied current pulse

Ie2= 3e-9;          % magnitude of applied current pulse

%% Parameters for the lif neurons
Vth = -0.054;               % threshold potential (to produce spike)
Vreset = -0.080;            % reset potential (post-spike)
Cm = 1e-9;               % total membrane capacitance
R = 10e6; % membrane resistance
E_L = -70e-3;               % leak potential 

E12_rev = -70e-3;
E21_rev = -70e-3;
G12 = 1e-6; 
G21 = 1e-6;
tau_syn = 10e-3; %time constant
tau_D = 1; %time constant for depression variable

I01 = 2e-9; %baseline current for cell 1
I02 = 2e-9; %baseline current for cell 2
p_r = 1;

% define white noise
sigm = 5e-12; %noise level
       
 
       
fliptimes= []; %record state switch times 
a=1; %keeps count of flips

% Set up applied current to cell 1
I1app=zeros(size(t));   
for i = 1:length(I1app)
    I1app(i) = I01;     % baseline (negative here)
end


% Set up applied current to cell 2
I2app=zeros(size(t));   
for i = 1:length(I2app)
    I2app(i) = I02;     % baseline (negative here)
end
    

%% Begin with all variables for cell 1
% Define the arrays
V1=zeros(size(t));          % membrane potential
V1(1) = E_L;                % set the initial value of voltage     

%% repeat with all variables for cell 2
% Define the arrays
V2=zeros(size(t));          % membrane potential
V2(1) = E_L;                % set the initial value of voltage     


D1 = ones(size(t));   %depression variable for cell1
D2 = ones(size(t));   %depression variable for cell 2
syn1=zeros(size(t));    % synaptic gating variable from spikes in cell 1
syn2=zeros(size(t));    % synaptic gating variable from spikes in cell 2
spikes1 = zeros(size(t));   % store spikes of cell 1
spikes2 = zeros(size(t));   % store spikes of cell 2
spike1now = 0;              % set to 1 when cell 1 is in a spike
spike2now = 0;              % set to 1 when cell 2 is in a spike

% nspikewidth is number of time-points to update conductance following a spike
nspikewidth = round(tau_syn*10/dt);  

 spikestate=0; %variable to keep track of which cell is firing. 0 if cell2 is firing
%% Now simulate trials, each with a different applied current

for i = 2:length(t)            % loop through all time points
   
    % next line: Forward Euler method to update membrane potential for
    % cell 1 first
    noise = randn(1)* sigm/sqrt(dt); %noise term to add to the total membrane potential. sigma * n(t)
 
    V1(i) = V1(i-1) + dt*((E_L - V1(i-1))/R + I1app(i-1) +G21*syn2(i-1)*(E21_rev -V1(i-1)) + noise)/Cm;
    %update gating variables and depression variable
    syn1(i) = syn1(i-1) + dt*(-syn1(i-1)/tau_syn);
    D1(i) = D1(i-1) + dt*((1-D1(i-1))/tau_D);

    if (V1(i) > Vth )            % if potential is above threshold and a spike is not recorded
       if(spikestate==0)
           spikestate = 1;
           fliptimes(a) = i;
           a = a+1;
       end
           spike1now = 1;                          % detect this spike
        spikes1(i) = 1;                         % record this spike time
        
        V1(i) = Vreset;          % reset the potential

        syn1(i) = syn1(i) + p_r*D1(i)*(1-syn1(i));
        D1(i) = D1(i)*(1-p_r);
        
    end

    if ( V1(i) < Vth - 0.010 )               % once the spike is over
        spike1now = 0;                          % set this to zero so we are ready for the next spike.
    end
    
    %repeat for cell 2
     V2(i) = V2(i-1) + dt*((E_L - V2(i-1))/R + I2app(i-1) +G12*syn1(i-1)*(E12_rev -V2(i-1)) + noise)/Cm;
    %update gating variables and depression variable
    syn2(i) = syn2(i-1) + dt*(-syn2(i-1)/tau_syn);
    D2(i) = D2(i-1) + dt*((1-D2(i-1))/tau_D);

    if V2(i) > Vth              % if potential is above threshold
        if(spikestate==1)
           spikestate = 0;
           fliptimes(a) = i;
           a = a+1;
       end
       spike2now = 1;                          % detect this spike
        spikes2(i) = 1;                         % record this spike time

        V2(i) = Vreset;          % reset the potential
         syn2(i) = syn2(i) + p_r*D2(i)*(1-syn2(i));
        D2(i) = D2(i)*(1-p_r);
    end
    if ( V2(i) < Vth - 0.010 )               % once the spike is over
        spike2now = 0;                          % set this to zero so we are ready for the next spike.
    end

   
end                    
    statedurations = diff(fliptimes); %produce a vector whose odd entries correspond to durations of one state and whose even entries correspond to durations of the other.
   display(a + " switches observed.")

spikesin2 = [];
spikesin1 = [];
b = 1; 
d = 1;
for k = 1:length(statedurations)
    if(mod(k,2) == 0)
        spikesin2(b) = statedurations(k)*dt;
        b = b + 1;
    else
         spikesin1(d) = statedurations(k)*dt;
        d = d + 1;
    end
end



   %plot histogram of the duration in each state
   figure(2)
h2 = histogram(spikesin2);
hold on
title("Histogram of the duration of spikes in cell 2")
xlabel('Duration (s)')
ylabel('count')
hold off

figure(3)
h1 = histogram(spikesin1);
hold on
title("Histogram of the duration of spikes in cell 1")
xlabel('Duration (s)')
ylabel('count')


    %% Finally plot membrane potentials as a function of time
figure(1)
% Membrane potential of cell 1 in the upper panel
subplot(2,2,1)      
plot(t,V1,'k', 'LineWidth',0.001)
xlabel('Time (sec)')
ylabel('Cell 1 Vm')
hold on
% Membrane potential of cell 2 in the lower panel
subplot(2,2,2)      
plot(t,V2,'k','LineWidth',0.001)
xlabel('Time (sec)')
ylabel('Cell 2 Vm')

%plot synaptic gating variables against time
subplot(2,2,3)      
plot(t,syn1,'k','LineWidth',0.001)
xlabel('Time (sec)')
ylabel('Cell 1 synaptic gating variable')

%plot synaptic gating variables against time
subplot(2,2,4)      
plot(t,syn2,'k','LineWidth',0.001)
xlabel('Time (sec)')
ylabel('Cell 2 synaptic gating variable')
hold off
