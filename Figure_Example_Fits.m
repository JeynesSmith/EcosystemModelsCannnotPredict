% This script creates a figure of the top 20% of fits to a specific 
% experimental data set defined by the parameter, DS. 

% clear all
clear all

% load datasets
load 'Experimental fitting results/Results' FittingResults
load 'Experimental data'/Analysis_Timeseries

% Select an exemplar data set
DS = 89;

% Call up the data set in question
TS = Analysis_Timeseries{DS,4};

% define an upper threshold of datasets to plot
Q_threshold = 0.2;

% delete the worst fits
Fits = squeeze(FittingResults(DS,:,:));
SSD = [Fits{:,2}]; Fits = Fits(1:length(SSD),:);
F = find(SSD > quantile(SSD,Q_threshold));
Fits(F,:) = [];

NumSpp = length(Fits{1,3}); % number of  species
NumSets = size(Fits,1); % number of parameter sets
disp(['Keeping below quantile ' num2str(Q_threshold)])
disp([num2str(length(Fits)) ' fits remain' ])

% Set up the figure
figure(3), clf; FS = 16; FSA = 12; MS = 15; LW = 2;
CL = get(gca,'colororder');

% simulate 10 random parameter sets
for p = 1:10
    % select a random parameter set
    R = randi(size(Fits,1));
    
    % Simulate the values for 25 timesteps
    T_obs = linspace(0,19.5,1000);
    
    % define network parameters
    Thisni = Fits{R,3}; % initial condition
    Thisr =  Fits{R,4}; % growth rate
    ThisA =  Fits{R,5}; % interaction matrix
    
    % simulate the network behaviour
    [t1,n1] = ode45(@(t,n)species_DE(t,n,ThisA,Thisr),T_obs,Thisni);
    
    % scale
    for i = 1:NumSpp
        N = max(TS(i,:));
        n1(:,i) = n1(:,i)./N;
    end
    
    % plot
    for s = 1:NumSpp
        subplot(NumSpp,1,s), hold on
        
        plot(t1,n1(:,s),'-' ,'color',CL(1,:)+0.1,'linewidth',1.5)
    end
end

% plot raw data 
for s = 1:NumSpp
    subplot(NumSpp,1,s), hold on
    
    plot([0:18],TS(s,:)./max(TS(s,:)),'o','markersize',5,'color','k','linewidth',1.5)
    ylim([0 2]); xlim([-0.25 max(T_obs)+0.5])
    ylabel(['$n_{' num2str(i) '}(t)$'],'fontsize',FS,'Interpreter','latex')
end
xlabel('Timesteps','fontsize',FS,'Interpreter','latex')

% save
% Make_TIFF('Figures/Figure_example_fits.tiff',[0 0 25 35])










