% Vincent Calia-Bogan
% Tutorial 5.3 -- bistability of two LIF neurons 
% % Spring 2024 
% Professor Paul Miller 
% Help Recieved: 
% Book Pages Refrenced: 
% Miller Codes Refrenced: 
% Functions called: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: Fig. 5.9a has a very similar codebase 
% decide if you want ISIs or not-- could be useful, might not be 

clear; 
% Base variables for the identical oscillators 
Cm = 1e-9; % membrane capacitance in (F) 
Rm = 10e6; % membrane resistance in (Ohm) 
E = -70e-3; % leak / reversal potential in (V) 
Vth = -54e-3; % threshold potential in (V) 
Vr = -80e-3; % reset potential in (V) 
% making sure synapses are identical and inhibitory
Erev12 = -70e-3; % reversal potential for one on two in (V) 
Erev21 = Erev12; % reversal potential for two on one in (V) 
G12 = 1e-6; % conductance for one on two in (S) 
G21 = G12; % conductance for two on one in (S) 
Tausyn = 10e-3; % time constant in synapse in(Sec) 
Iappbl = 2e-9; % baseline applied current in (A)
D1 = 1; % Synaptic depression for neuron one is one during this simulation 
D2 = D1; % synamptic depression is the same for both 
Pr = 1; % Synaptic depression additive bit (?)
TauD = 1; % Time constant for depression during the simulation -- pt. B

% universal vectors
tmax = 6; % maximum time for the sim 
dt = 1e-6; % dt-- make this smaller if needed 
tvec = 0:dt:tmax; % time vector 
sigma = 0e-12; % set initially to 0 pA * s .^ 0.5 -- input noise 
stdevnoise = (sigma / sqrt(dt)); % setting up standard deviation for rand. num. 
randmod = randn(size(tvec)) + 0 ; % mean of 0 
Nt = zeros(size(tvec)); % vector for Nt, the white noise being added (see tutorial 3.2)

% vectors for neuron 1 
Vvec1 = zeros(size(tvec)); % voltage for neuron 1
Vvec1(1) = E; % set initial membrane potential Eleak
Svec1 = zeros(size(tvec)); % synaptic gating variable for neuron 1 
Dvec1 = ones(size(tvec)) * D1; % depression for neuron 1 
Iapp1 = ones(size(tvec)) * Iappbl; % applied current for neuron 1
fr1 = zeros(size(tvec));  % firing rate for neuron 1
% adding additional current 
Iadd1 = 3e-9 ; % additional applied current to be added
Iaddpts1 = 0.1 / dt; % number of time points for adding extra current; 100 ms at the front
Iapp1(1:min(1 + Iaddpts1 - 1, numel(Iapp1))) = Iapp1(1:min(1 + Iaddpts1 - 1, numel(Iapp1))) + Iadd1; 
% adding additional applied current at the start

% vectors for neuron two 
Vvec2 = zeros(size(tvec)); % voltage for neuron 2
Vvec2(1) = E; % set initial membrane voltage to Eleak
Svec2 = zeros(size(tvec)); % synaptic gating variable for neuron 2 
Dvec2 = ones(size(tvec)) * D2; % depression for neuron 2 
Iapp2 = ones(size(tvec)) * Iappbl; % applied current for neuron 2
fr2 = zeros(size(tvec));  %firing rate for neuron 2 
% adding additional current 
Iadd2 = 3e-9 ; % additional applied current to be added
Iaddpts2 = 0.1 / dt; % number of time points for adding extra current; 100 ms at the front
Istart = length(Iapp2) / 2 - (Iaddpts2 / 2); % position for adding current to second neuron (halfway through)
Iapp2(round(Istart):round(Istart + Iaddpts2 - 1, numel(Iapp2))) = ...
Iapp2(round(Istart):round(Istart + Iaddpts2 - 1, numel(Iapp2))) + Iadd2; 
% adding additional applied current at the start


%% As accurate to the book as I can make it [ no work ] 
% NOTE: This write-up is dead accurate to the book so far as I can figure, but doesn't show anything 
% update: near as I can tell, this is identical to the stated equations in
% the book. Why this doesn't work is a mystery to me atm 
for i = 2:length(tvec)
    Nt(i) = (stdevnoise * (randmod(i))) + 0; % noisey input with correct standard deviation
    dV1Dt = ((E - Vvec1(i-1)) / Rm + G21 * Svec2(i-1) * (Erev21 - Vvec1(i-1)) + Iapp1(i-1) + Nt(i-1))/ Cm;
    Vvec1(i) = Vvec1(i-1) + dV1Dt * dt; 
    dS1Dt = - Svec1(i-1) / Tausyn;
    Svec1(i) = Svec1(i-1) + dS1Dt * dt; 
    %dD1Dt = (1 - Dvec1(i-1)) / TauD; later on
    %Dvec1(i) = Dvec1(i-1) + dD1Dt * dt; later on 
    if Vvec1(i) > Vth
        Vvec1(i) = Vr; 
        Svec1(i) = Svec1(i) + Pr * Dvec1(i) * (1-Svec1(i)); 
        % Dvec1(i) = Dvec1(i) * (1-Pr); later 
        fr1(i) = 1; % write a 1 to firing rate 
    end
    dV2Dt = ((E - Vvec2(i-1)) / Rm + G12 * Svec1(i-1) * (Erev12 - Vvec2(i-1)) + Iapp2(i-1) + Nt(i-1))/ Cm;
    Vvec2(i) = Vvec2(i-1) + dV2Dt * dt; 
    dS2Dt = - Svec2(i-1) / Tausyn;
    Svec2(i) = Svec2(i-1) + dS2Dt * dt; 
    %dD2Dt = (1 - Dvec2(i-1)) / TauD; later on
    %Dvec2(i) = Dvec2(i-1) + dD2Dt * dt; later on 
    if Vvec2(i) > Vth
        Vvec2(i) = Vr; 
        Svec2(i) = Svec2(i) + Pr * Dvec2(i) * (1-Svec2(i)); 
        % Dvec2(i) = Dvec2(i) * (1-Pr); later 
        fr2(i) = 1; % write a 1 to firing rate 
    end

end


%% Trying a different methodology-- exponential method of solve for voltage [ no work ] 
% NOTE: This does show oscillations, but doesn't show switching between
% neuron 1 and 2-- perhaps for the same reason as above. 

%neuron 1 extra vectors
Il1= zeros(size(tvec));       % leak current
Isyn1 = zeros(size(tvec)); % synaptic current from neuron 2 into neuron 1 
Itot1 = zeros(size(tvec)); 
state1 = 0; % create an int that gets iterated for each state
%TODO: perhaps another way to do this? 
% extra vectors for neuron 2
Il2= zeros(size(tvec));       % leak current
Isyn2 = zeros(size(tvec)); % synaptic current from neuron 2 into neuron 1 
Itot2 = zeros(size(tvec)); 
state2 = 0; % create an int that gets iterated for neuron 2 
% trying something 
nspikewidth = round(Tausyn/dt);  

for i = 2:length(tvec) % iterate through time
    Nt(i) = (stdevnoise * (randmod(i))) + 0; % noisey input with correct standard deviation
    % neuron 1
    Dvec1(i-1) = D1; % for part A; ensuring depression is 1 
    Il1(i) = (E-Vvec1(i-1))/Rm;  % Leak current from neuron 1
    Isyn1(i) = G21 * Svec2(i-1)*(E-Vvec1(i-1)); % synaptic current
    Itot1(i) = Il1(i-1)+ Isyn1(i-1) + Iapp1(i-1) + Nt(i-1); % sum of all current + additional noise
    Gtot1 = inv(Rm) + G21 * Isyn1(i-1); % total conductance for the cell 
    Vinf = ((E/Rm + G21 * Isyn1(i-1) * Erev21 + Iapp1(i-1))/Gtot1); % Vinf is steady state
    Vvec1(i) = Vinf - (Vinf-Vvec1(i-1))*exp(-dt*Gtot1/Cm);  % iterate Vmem via exponential method
    dS1Dt = - Svec1(i-1) / Tausyn ; % iterating synaptic gating variable 
    Svec1(i) = Svec1(i-1) + dS1Dt * dt; % iterating svec 
    % dD1Dt = (1-Dvec1(i-1)) / TauD ; % for part b only 
    % Dvec1(i) = Dvec1(i-1) + dD1Dt * dt ; % iterating depress, part b only
    if (Vvec1(i) > Vth ) && (state1 == 0)
        Vvec1(i) = Vr; % reset membrane voltage
        state1 = 1; % set state for when spike happens 
        fr1(i) = fr1(i) + 1; % record the spike at that time point-- better to just make this 1? 
      % Dvec1(i) = Dvec1(i) * (1-Pr); % for part B only
        Svec1(i) = Svec1(i) + Pr * D1 * (1-Svec1(i)); %issue: Both Svec stay at 0 ... 
    end
    if (Vvec1(i) < Vth - 0.010)
        state1 = 0; % setting state back to 0 when spike is done
    end
    % neuron 2
    Dvec2(i-1) = D2; % for part A; ensuring depression is 1 
    Il2(i) = (E-Vvec2(i-1))/Rm;  % Leak current from neuron 1
    Isyn2(i) = G12*Svec1(i-1)*(E-Vvec2(i-1)); % synaptic current
    Itot2(i) = Il2(i-1)+ Isyn2(i-1) + Iapp2(i-1) + Nt(i-1); % sum of all current + additional noise
    Gtot2 = inv(Rm)+ G12 * Isyn2(i-1); % total conductance for the cell 
    Vinf = (E/Rm + G12* Isyn2(i-1)*Erev12 + Iapp2(i-1))/Gtot2; % Vinf is steady state
    Vvec2(i) = Vinf - (Vinf-Vvec2(i-1))*exp(-dt*Gtot2/Cm);  % iterate Vmem via exponential method
    dS2Dt = - Svec2(i-1) / Tausyn ; % iterating synaptic gating variable 
    Svec2(i) = Svec2(i-1) + dS2Dt * dt; % iterating svec 
    % dD2Dt = (1-Dvec2(i-1)) / TauD ; % for part b only 
    % Dvec2(i) = Dvec2(i-1) + dD2Dt * dt ; % iterating depress, part b only
    if (Vvec2(i) > Vth ) && (state2 == 0)
        Vvec2(i) = Vr; % reset membrane voltage
        state2 = 1; % set state for when spike happens 
        fr2(i) = fr2(i) + 1; % record the spike at that time point-- better to just make this 1? 
        %Dvec1(i) = Dvec1(i) * (1-Pr); % for part B only
        Svec2(i) = Svec2(i) + Pr * D2 * (1-Svec2(i));
    end 
    if (Vvec2(i) < Vth - 0.010)
        state2 = 0; % setting state back to 0 when spike is done
    end
end
%% Forward euler with chunking according to current-- not working [ no work ] 
% simulation, pt.A, effects of synaptic depression are removed: 
% NOTE: This has exactly the same issue as the very top section-- again
% likely due to the same issue 

% no hoppin clue why this no work 
for i=2:length(tvec)
    Nt(i) = (stdevnoise * (randmod(i))) + 0; % noisey input with correct standard deviation
    % neuron 1 segment
    Dvec1(i-1) = D1; % for part A; ensuring depression is 1 
    Iv1 = (E - Vvec1(i-1)/ Rm); % current due to membrane voltage 
    Ig1 = G21 * Svec2(i-1) * (Erev21 - Vvec1(i-1)); % current due to condutance and synaptic gating var.
    dV1Dt = (Iv1 + Ig1 + Iapp1(i-1)) / Cm; % forward euler iterating neuron 1 voltage 
    Vvec1(i) = Vvec1(i-1) + dV1Dt * dt + Nt(i); % iterating voltage 
    dS1Dt = - Svec1(i-1) / Tausyn ; % iterating synaptic gating variable 
    Svec1(i) = Svec1(i-1) + dS1Dt * dt; % iterating svec 
    % dD1Dt = (1-Dvec1(i-1)) / TauD ; % for part b only 
    % Dvec1(i) = Dvec1(i-1) + dD1Dt * dt ; % iterating depress, part b only
    if (Vvec1(i) > Vth)
        Vvec1(i) = Vr; % reset membrane voltage 
        Svec1(i) = Svec1(i) + Pr * Dvec1(i) * (1-Svec1(i));
        % Dvec1(i) = Dvec1(i) * (1-Pr); % for part B only
        fr1(i) = fr1(i) + 1; % count up the firing rate for neuron one 
    end 
    % neuron 2 segment: 
    Dvec2(i-1) = D2; % for part A; ensuring depression is 1 
    Iv2 = (E-Vvec2(i-1) / Rm); % current due to membrane voltage 
    Ig2 = G12 * Svec2(i-1) * (Erev12 - Vvec2(i-1)); % current due to condutance and synaptic gating var.
    dV2Dt =(Iv2 + Ig2 + Iapp2(i-1) + Nt(i)) / Cm; % forward euler iterating neuron 1 voltage 
    Vvec2(i) = Vvec2(i-1) + dV2Dt * dt; % iterating voltage 
    dS2Dt = -Svec2(i-1) / Tausyn ; % iterating synaptic gating variable 
    Svec2(i) = Svec2(i-1) + dS2Dt * dt; % iterating svec 
    % dD2Dt = (1-Dvec2(i-1)) / TauD ; % for part b only 
    % Dvec2(i) = Dvec2(i-1) + dD2Dt * dt ; % iterating depress, part b only
    if (Vvec2(i) > Vth)
        Vvec2(i) = Vr; % reset membrane voltage 
        Svec2(i) = Svec2(i) + Pr * Dvec2(i) * (1-Svec2(i));
        % Dvec1(i) = Dvec1(i) * (1-Pr); % for part B only
        fr2(i) = fr2(i) + 1; % count up the firing rate for neuron one 
    end 
end 

