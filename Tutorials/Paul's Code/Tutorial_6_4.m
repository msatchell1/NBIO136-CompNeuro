% Tutorial_6_4.m  
% Code for contrast-invariant orientation selectivity
%
%   Three types of network are simulated. 
%   If question_part = 'A' A feedforward excitatory network is used.
%   If question_part = 'B' then feedforward cross-inhibition is added.
%   If question_part = 'C' then a recurrent network (both E and I) is used.
%
% This code is a solution to Tutorial 6.4 of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;                          % Clear all variables

question_part = 'C';            % Can be 'A' or 'B' or 'C'

dt=0.0001;                      % Time step in sec (0.1ms)
tmax =0.3;                      % Time step in sec (300ms)
t=0:dt:tmax;                    % Create time vector

Ncells = 50;                    % Number of cells around the ring

rE=zeros(length(t),Ncells);     % Rate matrix for all timesteps and all cells
rI=zeros(length(t),Ncells);     % Rate matrix for all timesteps and all cells
hE=zeros(1,Ncells);             % Applied thalamic input to all cells
IstimE=zeros(1,Ncells);         % Total current to all cells
IstimI=zeros(1,Ncells);         % Total current to all cells

%% Cortical network parameters
switch question_part
    case 'A'
        AE = 40;            % Maximum LGN input to E cells
        AI = 40;            % Maximum LGN input to I cells
        I0E = -10;          % Background input to E cells minus threshold
        I0I = -10;          % Background input to I cells minus threshold
        tauE = 0.010;       % Time constant for E cells
        tauI = 0.010;       % Time constant for I cells
        WEE0 = 0;           % Mean E to E connection weight
        WEI0 = 0;           % Mean E to I connection weight
        WIE0 = 0;           % Mean I to E connection weight
        WIEshift = pi;      % Difference in tuning preference for I cells connecting to E cells      
    case 'B'
        AE = 40;            % Maximum LGN input to E cells
        AI = 40;            % Maximum LGN input to I cells
        I0E = -5;           % Background input to E cells minus threshold
        I0I = -5;           % Background input to I cells minus threshold
        tauE = 0.010;       % Time constant for E cells
        tauI = 0.010;       % Time constant for I cells
        WEE0 = 0;           % Mean E to E connection weight
        WEI0 = 0;           % Mean E to I connection weight
        WIE0 = -1;          % Mean I to E connection weight
        WIEshift = pi;      % Difference in tuning preference for I cells connecting to E cells      
    case 'C'
        AE = 100;           % Maximum LGN input to E cells
        AI = 0;             % Maximum LGN input to I cells
        I0E = 2;            % Background input to E cells minus threshold
        I0I = 0.5;            % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 5;           % Mean E to E connection weight
        WEI0 = 3;           % Mean E to I connection weight
        WIE0 = -4;          % Mean I to E connection weight
        WIEshift = 0;       % Difference in tuning preference for I cells connecting to E cells
end

% Now produce all the within-network connections as cosine functions. Note
% that 1 is added to the cosine so the variation of the term within
% parentheses has a minimum of 0 and a maximum of 2. 
% The initial terms WEE0, WEI0, and WIE0 correspond then to the mean weight
% of that type of connection.
% The I-to-E connection has an extra "WIEshift" term which can be set to pi
% so that inhibition is from cells with opposite LGN input.
for cell1 = 1:Ncells         
    for cell2 = 1:Ncells            
        WEE(cell1,cell2) = WEE0*(1+cos(2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WEI(cell1,cell2) = WEI0*(1+cos(2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WIE(cell1,cell2) = WIE0*(1+cos(WIEshift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
    end
end

%% Stimulus input parameters
epsilon = 0.5;              % Variation of input across cells from 1-epsilon to 1+epsilon
cuecell = Ncells/2;         % Cell with peak of the LGN input
hE = zeros(1,Ncells);       % Vector for inputs to each E cell when input is on
hI = zeros(1,Ncells);       % Vector for inputs to each I cell when input is on
hstart = 0.0;               % Time to begin input
hend = tmax;                % Time to end input
contrasts = [0 0.25 0.5 0.75 1];    % Range of contrasts to use

%% Set up the plotting parameters 
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

% Now loop through the set of different contrasts
for trial = 1:length(contrasts) 
    c = contrasts(trial);               % Contrast to be used
    IstimE = zeros(1,Ncells);           % Define vector for input current to E cells
    IstimI = zeros(1,Ncells);           % Define vector for input current to E cells
    
    % Now set the input current that varies across cells
    for cell = 1:Ncells         
        hE(cell) = AE*c*(1 + epsilon*cos(2*pi*(cell-cuecell)/Ncells));
        hI(cell) = AI*c*(1 + epsilon*cos(2*pi*(cell-cuecell)/Ncells));
    end
    
    % Begin the time integration
    for i = 2:length(t)         
              
        % Only set the input current to each unit as the input current when
        % the stimulus is on from hstart to hend.
        if ( t(i) >= hstart ) && (t(i) <= hend )
            IstimE = hE;                        % Set input to E cells
            IstimI = hI;                        % Set input to I cells
        else
            IstimE = zeros(1,Ncells);           % Otherwise no input current
            IstimI = zeros(1,Ncells);           % Otherwise no input current
        end
        
        % Update rates of all excitatory units based on rates of all units
        % at the previous time point
        rE(i,:) = rE(i-1,:)*(1-dt/tauE) + ...
            dt*(IstimE + rE(i-1,:)*WEE + rI(i-1,:)*WIE + I0E)/tauE;
        
        % Update rates of all inhibitory units based on rates of all units
        % at the previous time point
        rI(i,:) = rI(i-1,:)*(1-dt/tauI) + ...
            dt*(IstimI + rE(i-1,:)*WEI +I0I)/tauI;

        rE(i,:) = max(rE(i,:),0);               % Rates cannot be less than 0
        rI(i,:) = max(rI(i,:),0);               % Rates cannot be less than 0
        
    end
    
    %% Now plot the results for the contrast used
    figure(1)
    if ( trial == 1 )
        clf;                    % Clear figure on first trial
    end
    
    subplot(2,1,1)
    % plot rate of all E cells at end of simulation
    plot(rE(length(t),:),'g')                
    hold on
    xlabel('cell index')
    ylabel('rate of E-neurons')
    
    subplot(2,1,2)    
    % plot rate of all I cells at end of simulation
    plot(rI(length(t),:),'g')                
    xlabel('cell index')
    ylabel('rate of I-neurons')
    
    figure(2)
    if ( trial == 1 )
        clf;                    % Clear figure on first trial
    end
    subplot(2,1,1)
    % plot rate of excitatory cuecell as a function of time
    plot(t,rE(:,cuecell))                    
    hold on
    % plot rate of E cell with null direction input as a function of time
    plot(t,rE(:,1+mod(cuecell+Ncells/2-1,Ncells)))                    
    xlabel('time')
    ylabel('rate of excitatory cell')
    
    
    subplot(2,1,2)
    % plot rate of inhibitory cuecell as a function of time
    plot(t,rI(:,cuecell))                    
    hold on
    % plot rate of I cell with null direction input as a function of time
    plot(t,rI(:,1+mod(cuecell+Ncells/2-1,Ncells)))                    
    xlabel('time')
    ylabel('rate of inhibitory cell')
    
    figure(3)
    if ( trial == 1 )
        clf;                    % Clear figure on first trial
    end
    
    % Plot rate of all excitatory cells normalized by the population mean
    plot(rE(end,:)/mean(rE(end,:)))
    hold on
end % Next contrast

figure(23)
mesh(WEE)
colormap(gray)
axis([0 50 0 50 0 0.25])
xlabel('Neuron i')
ylabel('Neuron j')
zlabel('WEE_{ij}')
