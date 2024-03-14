%tutorial 4.1 Q1c
% Sayaka (Saya) Minegishi
% 2/10/2024

clear
close all

dt = 1e-6;          % time-step for integration (sec)
tmax=0.35;           % maximum time of simulation (sec)
t=0:dt:tmax;            % vector of time points

%% Neuron parameters
E_L = -0.060;       % leak reversal potential (V)
E_Na = 0.045;       % reversal for sodium channels (V)
E_K = -0.082;       % reversal for potassium channels (V)
V0 = -0.0702; %membrane potential stabilizes at this voltage                 (V)

G_L = 30e-9;        % specific leak conductance (S)
G_Na = 12e-6;       % specific sodium conductance (S)
G_K = 3.6e-6;         % specific potassium conductance (S)

Cm = 100e-12;       % specific membrane capacitance (F)

Ntrials = 10; %number of current pulses
amplitude_iapp = 0.22e-9; %amplitude of applied current (A)
Ibase = 0;          % Baseline current before pulse
Iapp=Ibase*ones(size(t));   % Applied current has a baseline


min_delay_pulse = 5e-3; %minimum delay between pulses (s)
max_delay_pulse = 25e-3; %max delay between pulses (s)
deltat = (max_delay_pulse-min_delay_pulse)/Ntrials;
delay_btwn_pulse = [min_delay_pulse : deltat : max_delay_pulse]; %delay between pulses (s)

istart = 0.1;      % time applied current starts (sec)
ilength=5e-3;        % length of applied current pulse (sec)


V=zeros(size(t));           % membrane potential vector
V(1) = E_L;             % set the inititial value of voltage
    
n=zeros(size(t));       % n: potassium activation gating variable
n(1) = 0.35;            % start off near steady state when V is E_L
m=zeros(size(t));       % m: sodium activation gating variable
m(1) = 0.05;            % start off near steady state when V is E_L
h=zeros(size(t));       % h: sodim inactivation gating variable
h(1) = 0.75;            % start off near steady state when V is E_L

Itot=zeros(size(t));    % record the total current
I_Na=zeros(size(t));    % record sodium curret
I_K=zeros(size(t));     % record potassium current
I_L=zeros(size(t));     % record leak current

%% Now loop through trials
for trial = 1:Ntrials
    istart = istart + delay_btwn_pulse(trial); %update start time of pulse
    Ie= amplitude_iapp;           % New applied current each trial
    
    %% current clamp initialization
    
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        Iapp(i) = Ie;
    end
    
    
    for i = round(istart/dt)+1:round((istart+ilength)/dt) % now see how things change through time
        
        I_L(i) = G_L*(E_L-V(i));      % calculate leak current
        
        Vm = -70-1000*V(i);             % convert V to Hodgkin-Huxley form               % Vm is instantaneous voltage in mV
        
      %% Update all of the voltage-dependent rate constants for gating variables
        if ( Vm == -25 )
            alpha_m = 0.1/0.1;              % sodium activation rate
        else
            alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);   % sodium activation rate
        end
        beta_m = 4*exp(Vm/18);              % sodium deactivation rate
        
        alpha_h = 0.07*exp(Vm/20);          % sodium inactivation rate
        beta_h = 1/(1+exp((Vm+30)/10));     % sodium deinactivation rate
        
        if ( Vm == -10)
            alpha_n = 0.01/0.1;             % potassium activation rate
        else
            alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);  % potassium activation rate
        end
        beta_n = 0.125*exp((Vm)/80);        % potassium deactivation rate
        
        %% Use rate constants to evaluate time constants and instantaneous steady state
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        tau_m = 1e-3/(alpha_m+beta_m);         % sodium activation variable
        m_inf = alpha_m/(alpha_m+beta_m);   % sodium activation variable
        
        tau_h = 1e-3/(alpha_h+beta_h);         % sodium inactivation variable
        h_inf = alpha_h/(alpha_h+beta_h);   % sodium inactivation variable
        
        tau_n = 1e-3/(alpha_n+beta_n);         % potassium activation variable
        n_inf = alpha_n/(alpha_n+beta_n);   % potassium activation variable
        
        
        % Now update gating variables using the Forward Euler method
        if ( i > 1 )
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        end
        
        % Now update currents and membrane potential using Forward Euler
        % method
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % sodium current
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i));    % potassium current
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i);        % total current is sum of leak + active channels + applied current
        V(i+1) = V(i) + Itot(i)*dt/Cm;                  % Update the membrane potential, V.
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(1) % Figure 4.6 in the textbook
  
    
    %% Plot all applied currents on one subplot at the top
    subplot('Position',[0.15 0.86 0.8 0.13])
    plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
    ylabel('I_{app} (nA)')
    hold on
    
    %% Plot membrane potential traces one below the other by trial
    subplot('Position',[0.15 0.86-trial*0.16 0.8 0.12])
    plot(t(10:10:end),1000*V(10:10:end),'k');
    axis([0 tmax -85 40])
    
    % The command strcat allows concatenation of different pieces of text
    % (strings) to be used in the legend. The command num2str converts the
    % numerical value of Ie used in the trial into a string that can be
    % added to the legend
    legstring = strcat('Iapp =  ',num2str(Ie*1e9),' nA')
    legend(legstring);                  % add legend to the figure
    
    if ( trial == Ntrials )             % Only label the time axis once
        xlabel('Time, sec')
    end
    
    ylabel('V_{m} (mV)')
    set(gca,'YTick',[-80:20:40])
    set(gca,'YTickLabel',{'-80' '' '-40' '' '0' '' '20'})

   
end

%% The annotation command can add text to any location on the figure.
%  Here it is used for figure labels A-F.
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.79 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.63 0.05 0.05],'String','C','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.47 0.05 0.05],'String','D','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.31 0.05 0.05],'String','E','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.15 0.05 0.05],'String','F','LineStyle','none','FontSize',16,'FontWeight','Bold')

