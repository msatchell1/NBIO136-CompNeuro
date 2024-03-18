%tutorial 5.3 part a
%Sayaka (Saya) Minegishi
% minegishis@brandeis.edu
% Mar 12 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%% Simulation parameters (identical across cells)
dt = 0.0001;                % time-step
tmax = 6; %maximum simulation time (s)
t = 0:dt:tmax;               % vector of time-points (s)


istart1 = 0;          % time applied current starts (sec) for cell1
istart2 = tmax/2;          % time applied current starts (sec) for cell 2
ilength=0.1;    % length of applied current pulse (sec)

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
sigm = 0; %noise level
       
 
% Set up applied current to cell 1
I1app=zeros(size(t));   
for i = 1:round(istart1/dt)
    I1app(i) = I01;     % baseline (negative here)
end
for i=round(istart1/dt)+1:round((istart1+ilength)/dt) % make non-zero for duration of current pulse
    I1app(i) = Ie1 + I01;     % applied current (in fact zero)
end
for i = round((istart1+ilength)/dt):length(I1app)
    I1app(i) = I01; %return to baseline current
end

% Set up applied current to cell 2
I2app=zeros(size(t));   % Applied current vector
for i = 1:round(istart2/dt)
    I2app(i) = I02;
end
for i=round(istart2/dt)+1:round((istart2+ilength)/dt) % make non-zero for duration of current pulse
    I2app(i) = Ie2 + I02;
end
for i = round((istart2+ilength)/dt):length(I2app)
    I2app(i) = I02;
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


%% Now simulate trials, each with a different applied current

for i = 2:length(t)            % loop through all time points
    % next line: Forward Euler method to update membrane potential for
    % cell 1 first
    noise = randn(1)* sigm/sqrt(dt); %noise term to add to the total membrane potential. sigma * n(t)
 
    V1(i) = V1(i-1) + dt*((E_L - V1(i-1))/R + I1app(i-1) +G21*syn2(i-1)*(E21_rev -V1(i-1)) + noise)/Cm;
    %update gating variables and depression variable
    syn1(i) = syn1(i-1) + dt*(-syn1(i-1)/tau_syn);
    D1(i) = D1(i-1) + dt*((1-D1(i-1))/tau_D);

    if V1(i) > Vth              % if potential is above threshold
        spikes1(i) = 1;          % record the spike at that time-point
        V1(i) = Vreset;          % reset the potential

        syn1(i) = syn1(i) + p_r*D1(i)*(1-syn1(i));
        D1(i) = D1(i)*(1-p_r);
        
    end
    
    %repeat for cell 2
     V2(i) = V2(i-1) + dt*((E_L - V2(i-1))/R + I2app(i-1) +G12*syn1(i-1)*(E12_rev -V2(i-1)) + noise)/Cm;
    %update gating variables and depression variable
    syn2(i) = syn2(i-1) + dt*(-syn2(i-1)/tau_syn);
    D2(i) = D2(i-1) + dt*((1-D2(i-1))/tau_D);

    if V2(i) > Vth              % if potential is above threshold
        spikes2(i) = 1;          % record the spike at that time-point
        V2(i) = Vreset;          % reset the potential
         syn2(i) = syn2(i) + p_r*D2(i)*(1-syn2(i));
        D2(i) = D2(i)*(1-p_r);
    end

   
end                    
    
    %% Finally plot membrane potentials as a function of time
figure(1)
% Membrane potential of cell 1 in the upper panel
subplot(2,2,1)      
plot(t,V1,'k')
xlabel('Time (sec)')
ylabel('Cell 1 Vm')
hold on
% Membrane potential of cell 2 in the lower panel
subplot(2,2,2)      
plot(t,V2,'k')
xlabel('Time (sec)')
ylabel('Cell 2 Vm')

%plot synaptic gating variables against time
subplot(2,2,3)      
plot(t,syn1,'k')
xlabel('Time (sec)')
ylabel('Cell 1 synaptic gating variable')

%plot synaptic gating variables against time
subplot(2,2,4)      
plot(t,syn2,'k')
xlabel('Time (sec)')
ylabel('Cell 2 synaptic gating variable')
hold off
