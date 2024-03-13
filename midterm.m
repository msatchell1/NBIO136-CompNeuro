clear
dt = 0.0001; % time-step (0.1ms)
t = 0:dt:1; % 1 sec of time-points
GL = 20e-9; % Leak conductance (20nS)
EL = -0.065; % Leak potential (-65mV)
Vth = -0.048; % Threshold (-48mV)
Cm = 100e-12; % Capacitance (100pF)
Vreset = -0.070; % Reset potential after a spike (-70mV)
V = EL*ones(size(t)); % Membrane Potential
spikes = zeros(size(t));% Initialize spikes to zero at all times
Iapp = 350e-12*ones(size(t)); % 350pA applied current
Iapp(1:round(end/2)) = 0; % Applied current switches on at 0.5sec
for i = 1:length(t) % Simulate through time
    dVdt = ( (EL-V(i-1))*GL + Iapp(i) )/Cm; % dV/dt equation
    V(i) = V(i-1) + dVdt; % Update voltage
    if ( V(i) > Vth ) % When voltage is above threshold
        spikes(i) = 1; % record a spike
        V(i+1) = Vreset; % Reset membrane potential
    end
end
rate = sum(spikes) % No semi-colon so the rate is printed out