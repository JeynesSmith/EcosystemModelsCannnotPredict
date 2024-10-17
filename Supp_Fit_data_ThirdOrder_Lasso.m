% This script generates supplementary fits to microcosm data using the
% third-order interactions. This script can also be applied with lasso
% methods and a cutoff. This script only fits to unique combinations of
% species in the microcosm dataset.

LassoVals = [10, 1, 0.1, 0.01 0.001 0.0001 0 0 0 0];
CutOffVals = [0.1 0.01 0.001 0.0001 0 0 0 0 0 0];
for j = [10]%3
for i = [10]%2

Target = 1000; % number of simulations for each experiment
load 'Experimental data'\ProposedTargets_Unique TargetExperiments
FittingResults = cell(length(TargetExperiments),Target,8); % initialise storage
LassoParameter = LassoVals(i); % This parameter affects the strength of the lasso method which minimises paras
CutOffParameter = CutOffVals(j); % This parameter determines a cutoff where parameters become 0
save(['Experimental fitting results\Results_ThirdOrder_Lasso' num2str(i) '_Cuttoff' num2str(j)],'FittingResults','Target','TargetExperiments','LassoParameter','CutOffParameter', '-v7.3') % save parameters and initial storage

% loop over microcosm experiments
parfor Index = 1:length(TargetExperiments)
    disp(Index)
    [A] = sub_Multiple_Experimental_Fits(TargetExperiments(Index),Target,LassoParameter,CutOffParameter); % run the fitting procedure
    FittingResults(Index,:,:) = A;
end
pause(1)
disp('Finished a round. Saving...')
save(['Experimental fitting results\Results_ThirdOrder_Lasso' num2str(i) '_Cuttoff' num2str(j)],'FittingResults','Target','TargetExperiments','LassoParameter','CutOffParameter', '-v7.3') % save parameters and initial storage
end
end

% fitting procedure
function [TempStorage] = sub_Multiple_Experimental_Fits(DS,Target,LassoParameter,CutOffParameter)

warning('off','all')

% Parameters
A_amp = 10; % This sets the range of the initial guesses. Larger values search more broadly

% Load experimental data
load 'Experimental data'\Analysis_Timeseries Analysis_Timeseries
TS = Analysis_Timeseries{DS,4}; % normalised Dataset
TS(TS==0)=1e-4;

% Let's define a "poor fit" (as opposed to a "good fit") as a model that
% does worse than a zero-order approximation (e.g., that n(t) = mean(n(t)))
Bar = ones(size(TS,1),size(TS,2)).*mean(TS,2); 
SSD_bar = sum(((TS - Bar)./TS).^2,'all') + sum(mean(TS,2));
if 1e6<SSD_bar
    disp('SSD Bar is too large')
    SSD_bar = 1e6;
end

NumSpp = size(TS,1); % How many species are there in this dataset
Tr = [1:size(TS,2)]'; % Create a time vector "Tr"
nr = TS'; % Transpose the abundance vector "nr"
n_max = 5.*[max(nr,[],1)]'; % upper bound for species abundances

% define the indices and matrix dimensions for higher order interactions
Indices = cell(NumSpp,NumSpp-2);
AllSpecies = 1:NumSpp;
% loop over species
for i = 1:NumSpp
    NotSpecies = AllSpecies; NotSpecies(i)=[];
    MatWidth = 0;
    for j = 2
        Indices{i,j-1} = nchoosek(NotSpecies,j)';
        MatWidth = MatWidth + size(Indices{i,j-1},2);
    end
end

% Define the upper and lower bounds of the search algorithm for all the
% values %%%%%%%%%%%%%%%%% must be the same for all
LB = -A_amp.*ones(NumSpp*(NumSpp+2 + MatWidth),1); % lower bound
UB = +A_amp.*ones(NumSpp*(NumSpp+2 + MatWidth),1); % upper bound

% No negative initial conditions OR growth rates allowed in the search
LB(1:2*NumSpp) = 0;
% LB(1:NumSpp) = 0; % allow negative growth

% Only negative diagonals allowed in the search
UB(2*NumSpp+1:NumSpp+1:(NumSpp*(NumSpp+2))) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check  
UB(1:NumSpp) = n_max; % upper bound for initial conditions

% Don't put outputs on the screen
options = optimset('Display','off');

% Search from random initial conditions until we have 100 decent fits to the data
count = 0; % initialise count
NumSets = 0; % initilise number of sets found
TempStorage = cell(Target,8); % allocate storage
while NumSets < Target % loop until Target is reached
    % update count, and display progress
    count = count + 1;
    
    % Randomly choose an initial parameter set, and solve the equations
    r_g = A_amp.*rand(NumSpp,1); % initial guesses for the growth rates
    % r_g = 2.*A_amp.*(rand(NumSpp,1)-0.5); % initial guesses for the growth rates % allow negative growth
    A_g = 2.*A_amp.*(rand(NumSpp)-0.5); % initial guesses for the interaction strengths
    A_g(1:NumSpp+1:end) = A_amp.*-1.*rand(NumSpp,1); % Make diagonal negative
    ni_g = abs(nr(1,:).*(randn([1,NumSpp])*0.5 + 1)); % initial conditions, based on observed
    B_g = 2.*A_amp.*(rand(NumSpp,MatWidth)-0.5); % initial guesses for the interaction strengths
    
    % Define the initial search location in parameter space as a 1D vector
    V0 = [ni_g(:); r_g(:); A_g(:); B_g(:)];
    
    % Run FMINCON with the inputs defined above
    [V_fit,SSD] = fmincon(@(x) FitNoisyData(x,Tr,nr,n_max,NumSpp,LassoParameter,Indices,MatWidth,CutOffParameter),V0,[],[],[],[],LB,UB,[],options);

    % We're only going to keep the best fit if it's better than the zero-order approximation
    if SSD < SSD_bar
        
        NumSets = NumSets + 1; % We've found one more best-fit dataset

        % % Add zeros
        % ZeroInd = find(abs(V_fit(NumSpp+1:end))<CutOffParameter);
        % if ~isempty(ZeroInd)
        %     V_fit(NumSpp+ZeroInd) = 0;
        % end
        
        % Store the details of the local minima
        TempStorage{NumSets,1} = V_fit; % Fitting vector
        TempStorage{NumSets,2} = SSD; % Fit statistic
        TempStorage{NumSets,3} = V_fit(1:NumSpp); % Initial conditions
        TempStorage{NumSets,4} = V_fit(NumSpp+1:2*NumSpp); % Growth rates
        TempStorage{NumSets,5} = reshape(V_fit(2*NumSpp+1:(NumSpp*(2+NumSpp))),NumSpp,NumSpp); % Interaction matrix
        TempStorage{NumSets,8} = reshape(V_fit((NumSpp*(2+NumSpp)+1):end),NumSpp,MatWidth); % Higher order interaction matrix
        % Simulate the model over the observation period
        opts = odeset('RelTol',1e-3,'Events',@(t,n) AbundEvent(t,n,n_max));
        [t_s,n_s] = ode45(@(t,n) species_DE(t,n,TempStorage{NumSets,5},TempStorage{NumSets,4},Indices,NumSpp,MatWidth,TempStorage{NumSets,8}),linspace(0,19,100),TempStorage{NumSets,3},opts);
        TempStorage{NumSets,6} = t_s; % Timevector
        TempStorage{NumSets,7} = n_s; % Best-fit modeled abundance vector
        disp(['For dataset #' num2str(DS) ': Number ' num2str(NumSets) ' found with ' num2str(sum(V_fit(NumSpp+1:end))==0) ' zeros and SSD=' num2str(SSD) '; ' num2str(count) ' searches.'])
    % else
    %     disp('Failed optimisation')
    end
end
end
% data fiting function
function [SSD,t_q,n_q] = FitNoisyData(V,t_n,n_n,n_max,NumSpp,parameter,Indices,MatWidth,CutOffParameter)

% % Break V vector into the contributing arrays
% ni_search = V(1:NumSpp); % initial condition
% r_search = V(NumSpp+1:2*NumSpp); % growth rates
% A_search = reshape(V(2*NumSpp+1:end),NumSpp,NumSpp); % interaction matrix

% % Add zeros
% ZeroInd = find(abs(V(NumSpp+1:end))<CutOffParameter);
% if ~isempty(ZeroInd)
%     V(NumSpp+ZeroInd) = 0;
% end

% solve GLV with parameters
opts = odeset('RelTol',1e-3,'Events',@(t,n) AbundEvent(t,n,n_max));
[t_q,n_q] = ode45(@(t,n) species_DE(t,n,reshape(V(2*NumSpp+1:(NumSpp*(2+NumSpp))),NumSpp,NumSpp),V(NumSpp+1:2*NumSpp),Indices,NumSpp,MatWidth,reshape(V((NumSpp*(2+NumSpp)+1):end),NumSpp,MatWidth)),t_n,V(1:NumSpp),opts);
n_q(n_q<0) = 0; % ensure all abundances are non-negative

% calculate SSD
if size(n_n,1) ~= size(n_q,1)  % if vector lengths do not match (because of parameter problems e.g. infinite growth/death), set SSD to be large
    SSD = 1e6;
else
    SSD = sum(((n_n-n_q)./n_n).^2,'all');% + parameter*sum(abs(V(NumSpp+1:end)));
%     SSD = sum(abs((n_n-n_q)./n_n),'all') + parameter*sum(abs(V(NumSpp+1:end)));
end
end

function [value, isterminal, direction] = AbundEvent(t,n,n_max)
value = any(n<0) || any((n_max - n)<0);
isterminal = 1;
direction = 0;
end

function dndt = species_DE(t,n,A,r,Ind,NS,MW,B)
BM = zeros(NS,MW);
for i = 1:NS
    temp = [];
    for j = 1:NS-2
        temp = [temp, prod(n(Ind{i,j}),1)];
    end
    BM(i,:) = temp;
end
% update ODE system based on the generalised Lotka-Volterra equations
dndt = n.*r + (A*n).*n + n.*sum(B.*BM,2);
end