% This script generates 500 fitted datasets to the 111 experimental microcosm 
% datasets. The fitted data is saved as FittingResults.mat. The number of
% fitted datasets are specified by Target, and the expereimental data which
% is fitted is denoted by DS. Parameters for a generalised Lotka-Volterra
% model are fitted to data using fmincon and rejected if they produce a
% sum-of-squared deviations greater than a null model (constant abundance =
% average abundance for every species).

Target = 500; % number of simulations for each experiment
FittingResults = cell(111,Target,7); % initialise storage
save 'Experimental fitting results'\Results FittingResults Target % save parameters and initial storage

% loop over microcosm experiments
for DS = 1:111
    disp(DS) 
    sub_Multiple_Experimental_Fits(DS) % run the fitting procedure
end

% fitting procedure
function sub_Multiple_Experimental_Fits(DS)

warning('off','all')

% Load experimental data
load 'Experimental data'\Analysis_Timeseries Analysis_Timeseries
load 'Experimental fitting results'\Results FittingResults Target
TS = Analysis_Timeseries{DS,4};

% Let's define a "poor fit" (as opposed to a "good fit") as a model that 
% does worse than a zero-order approximation (e.g., that n(t) = mean(n(t)))
MTS = TS'; MTS = MTS./repmat(mean(MTS),19,1);
SSD_bar = sum(sum((MTS - ones(size(MTS))).^2));
SSD_bar = SSD_bar;

NumSpp = size(TS,1); % How many species are there in this dataset
Tr = [1:size(TS,2)]'; % Create a time vector "Tr"
nr = TS'; % Transpose the abundance vector "nr"

% Search from random initial conditions until we have 100 decent fits to the data
count = 0; % initialise count
NumSets = 0; % initilise number of sets found
while NumSets < Target % loop until Target is reached
    % update count, and display progress
    count = count + 1; 
    disp(['For dataset #' num2str(DS) ':' num2str(count) ' searches; ' num2str(NumSets) ' found.'])
    
    % Randomly choose an initial parameter set, and solve the equations
    A_amp = 3; % This sets the range of the initial guesses. Larger values search more broadly
    r_g = A_amp.*rand(NumSpp,1); % initial guesses for the growth rates
    A_g = A_amp.*(rand(NumSpp)-0.5); % initial guesses for the interaction strengths
    A_g(1:NumSpp+1:end) = A_amp.*-1.*rand(NumSpp,1); % Make diagonal negative
    ni_g = abs(nr(1,:).*(randn([1,NumSpp])*0.5 + 1)); % initial conditions, based on observed
    
    % Define the initial search location in parameter space as a 1D vector
    V0 = [ni_g(:); r_g(:); A_g(:)];
    
    % Define the upper and lower bounds of the search algorithm for all the values
    BB = 50; % Magnitude of bounds
    LB = -BB.*ones(size(V0)); % lower bound
    UB = +BB.*ones(size(V0)); % upper bound
    
    % No negative initial conditions OR growth rates allowed in the search
    LB(1:2*NumSpp) = 1e-4;
    
    % Only negative diagonals allowed in the search
    UB(2*NumSpp+1:NumSpp+1:end) = -1e-3;
    
    % Don't put outputs on the screen
    options = optimset('Display','off');
    
    % Run FMINCON with the inputs defined above
    [V_fit,SSD] = fmincon(@(x) FitNoisyData(x,Tr,nr),V0,[],[],[],[],LB,UB,[],options);
    
    % We're only going to keep the best fit if it's better than the zero-order approximation
    if SSD < SSD_bar
        
        NumSets = NumSets + 1; % We've found one more best-fit dataset

        % Store the details of the local minima 
        FittingResults{DS,NumSets,1} = V_fit; % Fitting vector
        FittingResults{DS,NumSets,2} = SSD; % Fit statistic
        FittingResults{DS,NumSets,3} = V_fit(1:NumSpp); % Initial conditions
        FittingResults{DS,NumSets,4} = V_fit(NumSpp+1:2*NumSpp); % Growth rates
        FittingResults{DS,NumSets,5} = reshape(V_fit(2*NumSpp+1:end),NumSpp,NumSpp); % Interaction matrix

        % Simulate the model over the observation period
        [t_s,n_s] = ode45(@(t,n)species_DE(t,n,FittingResults{DS,NumSets,5},FittingResults{DS,NumSets,4}),linspace(0,19,100),FittingResults{DS,NumSets,3});
        FittingResults{DS,NumSets,6} = t_s; % Timevector
        FittingResults{DS,NumSets,7} = n_s; % Best-fit modeled abundance vector
    end
end
% save the experimental fit after the target has been reached
save 'Experimental fitting results'\Results FittingResults Target
end
% data fiting function
function [SSD,t_q,n_q] = FitNoisyData(V,t_n,n_n)

NumSpp = size(n_n,2); % number of species

% Break V vector into the contributing arrays
ni_search = V(1:NumSpp); % initial condition
r_search = V(NumSpp+1:2*NumSpp); % growth rates
A_search = reshape(V(2*NumSpp+1:end),NumSpp,NumSpp); % interaction matrix

% solve GLV with parameters
opts = odeset('RelTol',1e-3);
[t_q,n_q] = ode45(@(t,n)species_DE(t,n,A_search,r_search),t_n,ni_search,opts);
n_q(n_q<0) = 0; % ensure all abundances are non-negative

% calculate SSD
if length(n_n(:))-length(n_q(:)) ~= 0 % if vector lengths do not match (because of parameter problems e.g. infinite growth/death), set SSD to be large
    SSD = 1e5;
else
    SSD = sum(sum((n_n-n_q).^2));
end

if max(sum(n_q < 1e-5)) > 8 % all speies went extinct, then make SSD very large to avoid extinction models
    SSD = 1e5;
end
end















