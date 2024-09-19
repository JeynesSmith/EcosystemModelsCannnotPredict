% This script generates 500 fitted datasets to the 111 experimental microcosm
% datasets. The fitted data is saved as FittingResults.mat. The number of
% fitted datasets are specified by Target, and the expereimental data which
% is fitted is denoted by DS. Parameters for a generalised Lotka-Volterra
% model are fitted to data using fmincon and rejected if they produce a
% sum-of-squared deviations greater than a null model (constant abundance =
% average abundance for every species).

CutOffVals = [0.1 0.01 0.001 0.0001];

for i = 2:4

CutOffParameter = CutOffVals(i);
Target = 1000; % number of simulations for each experiment
load 'Experimental data'/ProposedTargets_Unique TargetExperiments
FittingResults = cell(length(TargetExperiments),Target,7); % initialise storage
save(['FitResults/Results_BAIC_Cuttoff' num2str(i) '_AllR'],'FittingResults','Target','TargetExperiments','CutOffParameter', '-v7.3') % save parameters and initial storage

% loop over microcosm experiments
parfor Index = 1:length(TargetExperiments)
    disp(Index)
    [A] = sub_Multiple_Experimental_Fits(TargetExperiments(Index),Target,CutOffParameter); % run the fitting procedure
    FittingResults(Index,:,:) = A;
end
pause(1)
disp('Finished a round. Saving...')
save(['FitResults/Results_BAIC_Cuttoff' num2str(i) '_AllR'],'FittingResults','Target','TargetExperiments','CutOffParameter', '-v7.3') % save parameters and initial storage
end

% fitting procedure
function [TempStorage] = sub_Multiple_Experimental_Fits(DS,Target,CutOffParameter)

warning('off','all')

% Parameters
A_amp = 10; % This sets the range of the initial guesses. Larger values search more broadly
SigmaMax = 10; % maximum standard deviation

% Load experimental data
load 'Experimental data'/Analysis_Timeseries Analysis_Timeseries
TS = Analysis_Timeseries{DS,4}; % normalised Dataset
TS(TS==0)=1e-4; %%%%%%%%%% do to all

% Let's define a "poor fit" (as opposed to a "good fit") as a model that
% does worse than a zero-order approximation (e.g., that n(t) = mean(n(t)))
Bar = ones(size(TS,1),size(TS,2)).*mean(TS,2); %%%%%%%%%%%%% do for all
SSD_bar = sum(((TS - Bar)./TS).^2,'all') + sum(mean(TS,2));
if 1e6<SSD_bar
    disp('SSD Bar is too large')
    SSD_bar = 1e6;
end

NumSpp = size(TS,1); % How many species are there in this dataset
Tr = [1:size(TS,2)]'; % Create a time vector "Tr"
nr = TS'; % Transpose the abundance vector "nr"
n_max = 5.*[max(nr,[],1)]'; % upper bound for species abundances

% Define the upper and lower bounds of the search algorithm for all the values
LB = -A_amp.*ones(NumSpp*(NumSpp+2)+1,1); % lower bound
UB = +A_amp.*ones(NumSpp*(NumSpp+2)+1,1); % upper bound

% No negative initial conditions OR growth rates allowed in the search
% LB(1:2*NumSpp) = 0;
LB(1:NumSpp) = 0; % allow negative growth
LB(end) = 1e-6; % make minimum SD very small

% Only negative diagonals allowed in the search
UB(2*NumSpp+1:NumSpp+1:end-1) = 0;
UB(1:NumSpp) = n_max; % upper bound for initial conditions
UB(end) = SigmaMax; % max standard deviation

% Don't put outputs on the screen
options = optimset('Display','off');

% Search from random initial conditions until we have 100 decent fits to the data
count = 0; % initialise count
NumSets = 0; % initilise number of sets found
TempStorage = cell(Target,7); % allocate storage
while NumSets < Target % loop until Target is reached
    % update count, and display progress
    count = count + 1;
    
    % Randomly choose an initial parameter set, and solve the equations
    % r_g = A_amp.*rand(NumSpp,1); % initial guesses for the growth rates
    r_g = 2.*A_amp.*(rand(NumSpp,1)-0.5); % initial guesses for the growth rates
    A_g = 2.*A_amp.*(rand(NumSpp)-0.5); % initial guesses for the interaction strengths
    A_g(1:NumSpp+1:end) = A_amp.*-1.*rand(NumSpp,1); % Make diagonal negative
    ni_g = abs(nr(1,:).*(randn([1,NumSpp])*0.5 + 1)); % initial conditions, based on observed
    sigma_g = rand*SigmaMax*4/5 + SigmaMax/5; % this ensures that SD isn't too small initially
    
    % Define the initial search location in parameter space as a 1D vector
    V0 = [ni_g(:); r_g(:); A_g(:); sigma_g];
    
    % Run FMINCON with the inputs defined above
    [V_fit,SSD] = fmincon(@(x) FitNoisyData(x,Tr,nr,n_max,NumSpp,CutOffParameter),V0,[],[],[],[],LB,UB,[],options);
    
    % We're only going to keep the best fit if it's better than the zero-order approximation
    if SSD < SSD_bar
        
        NumSets = NumSets + 1; % We've found one more best-fit dataset

        % Add zeros
        ZeroInd = find(abs(V_fit(NumSpp+1:end-1))<CutOffParameter);
        if ~isempty(ZeroInd)
            V_fit(NumSpp+ZeroInd) = 0;
        end
        
        % Store the details of the local minima
        TempStorage{NumSets,1} = V_fit; % Fitting vector
        TempStorage{NumSets,2} = SSD; % Fit statistic
        TempStorage{NumSets,3} = V_fit(1:NumSpp); % Initial conditions
        TempStorage{NumSets,4} = V_fit(NumSpp+1:2*NumSpp); % Growth rates
        TempStorage{NumSets,5} = reshape(V_fit(2*NumSpp+1:end-1),NumSpp,NumSpp); % Interaction matrix
        % Simulate the model over the observation period
        opts = odeset('RelTol',1e-3,'Events',@(t,n) AbundEvent(t,n,n_max));
        [t_s,n_s] = ode45(@(t,n) species_DE(t,n,TempStorage{NumSets,5},TempStorage{NumSets,4}),linspace(0,19,100),TempStorage{NumSets,3},opts);
        TempStorage{NumSets,6} = t_s; % Timevector
        TempStorage{NumSets,7} = n_s; % Best-fit modeled abundance vector
        disp(['For dataset #' num2str(DS) ': Number ' num2str(NumSets) ' found with ' num2str(sum(V_fit(NumSpp+1:end))==0) ' zeros and SSD=' num2str(SSD) '; ' num2str(count) ' searches.'])
    end
end
end
% data fiting function
function [SSD,t_q,n_q] = FitNoisyData(V,t_n,n_n,n_max,NumSpp,CutOffParameter)

% % Break V vector into the contributing arrays
% ni_search = V(1:NumSpp); % initial condition
% r_search = V(NumSpp+1:2*NumSpp); % growth rates
% A_search = reshape(V(2*NumSpp+1:end-1),NumSpp,NumSpp); % interaction matrix

% Add zeros
ZeroInd = find(abs(V(NumSpp+1:end-1))<CutOffParameter);
if ~isempty(ZeroInd)
    V(NumSpp+ZeroInd) = 0;
end

% solve GLV with parameters
opts = odeset('RelTol',1e-3,'Events',@(t,n) AbundEvent(t,n,n_max));
[t_q,n_q] = ode45(@(t,n) species_DE(t,n,reshape(V(2*NumSpp+1:end-1),NumSpp,NumSpp),V(NumSpp+1:2*NumSpp)),t_n,V(1:NumSpp),opts);
n_q(n_q<0) = 0; % ensure all abundances are non-negative

% calculate SSD
if size(n_n,1) ~= size(n_q,1)  % if vector lengths do not match (because of parameter problems e.g. infinite growth/death), set SSD to be large
    SSD = 1e6;
else
    SSD = 2*(NumSpp*(1+NumSpp)-length(ZeroInd))-2*sum(log(normpdf(n_q(:),n_n(:),V(end)))); % log likelihood
end
end

function [value, isterminal, direction] = AbundEvent(t,n,n_max)
value = any(n<0) || any((n_max - n)<0);
isterminal = 1;
direction = 0;
end

function dndt = species_DE(t,n,A,r,Press_amount,Press_species)
% this function defines the generalised Lotka-Volterra equations for a
% specific network. This includes a press perturbation applied to a species
%
% % t = time variable for ODE system
% % n = current abundances (vertical vector)
% % A = interaction matrix (square matrix with length = length(n))
% % r = growth rate vector (vertical vector with length = length(n))
% % Press_amount = strength of a press perturbation applied to a species
% % Press_species = species which the perturbation is applied to

% if no perturbation is defined, then apply a perturbation with zero
% strength to the first species (i.e. no effect)
if nargin == 4
    Press_amount = 0;
    Press_species = 1;
end

% update ODE system based on the generalised Lotka-Volterra equations
dndt = n.*r + (A*n).*n;

% Apply the perturbation proportional to species' abundance
dndt(Press_species) = dndt(Press_species) - Press_amount.*n(Press_species);
end