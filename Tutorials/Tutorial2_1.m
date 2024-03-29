% From "An Introductory Course in Computational Neuroscience"
% by Paul Miller
%
% Tutorial 2.1: The f-I curve of the leaky integrate-and-fire neuron



% Define parameters
El = -70; % mV
Rm = 5; % mega Ohms
Cm = 2; % nF
Vth = -50; % mV
Vreset = -65; % mV

% Create time vector
dt = 0.0001; % seconds
times = 0:dt:2;

% voltage vector of same size
V = zeros(size(times));

% Define initial condition of V(t=0) = El.
V(1) = El;

% Applied current values
currents = 0:1:500;
% resultiung firing rates
rates = zeros(size(currents));

for j = 1:numel(currents)
    % applied current vector
    I0 = currents(j);
    Iapp = zeros(size(times)) + I0;

    numSpikes = 0;
    
    for i = 2:numel(times)
        dV = ((El-V(i-1))/Rm + Iapp(i-1))*dt/Cm;
        V(i) = V(i-1) + dV;
    
        if V(i) > Vth % reset voltage if above threshold
            V(i) = Vreset;
            numSpikes = numSpikes + 1;
        end
    
    end
    
    rates(j) = numSpikes/(numel(times)*dt);

end


figure
plot(currents,rates)