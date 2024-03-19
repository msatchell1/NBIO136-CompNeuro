% Tutorial_9_1.m
%
% Generate multi-neuron data with only two degrees of freedom and analyze
% using PCA.
%
% Noise is added to the firing rate of each neuron, whose firing rates
% otherwise are produced by a linear combination of the two underlying
% degrees of freedom.
%
% The goal is to reconstruct the rates using the first two principle
% components alone and together.
%
% Code for the computational Tutorial 9.1 of the textbook
% "An Introductory Course in Computational Neuroscience"
% by Paul Miller, Brandeis University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
question_part = 1;      % Can be 1, 2, or 3 for different questions
%rng(3)
Ncells = 50;            % Number of distinct neural firing rates to use

dt = 0.001;             % dt in sec
tmax = 10;              % tmax in sec
t = [0:dt:tmax]';          % time vector
Nt = length(t);         % Number of time points

rate = zeros(Nt,Ncells); % firing rates of all cells is the main matrix

%% The following section generates dummy data from particular inputs
period = 2;
omega = 2*pi/period;

A0 = 50;
A1 = 20;
A2 = 10;
switch question_part
    case 1      
        %  In Question 1, two out of phase oscillations, which correspond 
        %  to "circular" or "ellipitical" motion, comprise the two inputs
        Input1 = A1*sin(omega*t);
        Input2 = A2*cos(omega*t);
        
    case 2
        %  In Question 2, two oscillations, one of double the frequency of 
        %  the first, comprise the two inputs.
        Input1 = A1*sin(2*omega*t);        
        Input2 = A2*cos(omega*t);
              
    case 3
        % In Question 3 a "ramp" has a brief "signal" superposed on it

        % Input 1 is a ramp of constant gradient
        Input1 = A1*t;

        % Input2 is a half-period of a sine wave within the "trial"
        ton1 = 2*period;
        toff1 = ton1+period/2;
        non1 = round(ton1/dt);
        noff1 = round(toff1/dt);

        Input2 = [zeros(non1,1); A2*sin(2*omega*t(non1+1:noff1)); zeros(Nt-noff1,1)];        
        
end

noise = 20;         % Level of random noise to be added to each cell


% The following weights make each neuron have a random linear combination 
% of content from each of the two inputs
W0 = 2+randn(1,Ncells);     % Constant, input-independent background rate
W1 = randn(1,Ncells);       % Component of Input 1 for each cell
W2 = randn(1,Ncells);       % Component of Input 2 for each cell


% The rate matrix has a row for each cell with time in columns
% This set of firing rates with added noise really only has two major
% signals in it, "Input1" and "Input2"
rate = A0*ones(Nt,1)*W0 + Input1*W1 + Input2*W2 + noise*randn(size(rate));

% Do PCA next. Each column should be a variable (i.e. a different cell) and
% each row an observation (here a time point). Need transpose of rate to do
% this correctly as rate stores the neurons as rows.
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(rate);

%% In this example the principal components are related to the weight matrices.
% Note that the sign (whether positively or negatively correlated) can not
% be determined from the principal component vectors in COEFF which could as
% easily point in the opposite directions.
figure(1)
plot(W1,COEFF(:,1),'x')

figure(2)
plot(W2,COEFF(:,2),'x')

% SCORE tells how much of each principal component is present at each point
% in time. If we simply multiplied COEFF*SCORE' we would recover the
% original data offset by the mean. By only taking the first two columns of COEFF (the first
% two principal components) and multiplying out the first two columns of
% SCORE (the first two rows of SCORE') we reduce the original 50-dimensions
% of data into a 2D data set.
COEFFT = COEFF';
newrates = SCORE(:,1:2)*COEFFT(1:2,:) + ones(Nt,1)*MU;

% Now we compare plots of the original rates and the 2D reduced rates.
figure(3)
subplot(2,1,1)
plot(rate(:,1),rate(:,2),'o')
subplot(2,1,2)
plot(newrates(:,1),newrates(:,2),'o')

figure(4)
subplot(2,1,1)
plot(t,rate(:,1))
subplot(2,1,2)
plot(t,newrates(:,2))


% Explained is % variance explained by successive components
figure(5)
plot(EXPLAINED,'--o')

figure(6)
subplot(2,1,1)
plot(t,SCORE(:,1))
subplot(2,1,2)
plot(t,SCORE(:,2))

