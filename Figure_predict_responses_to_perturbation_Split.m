function Figure_predict_responses_to_perturbation_Split()
% This function creates figures to examine the effect of perturbing a
% species within the network. We examine two comparisons: (i) whether a
% perturbation results in a better or worse outcome for a species compared
% to no perturbation, and (ii) whether a perturbation will increase or
% decrease a species' abundance.

% load data
load 'Experimental fitting results'/Results FittingResults

% open figures
figure(4), clf; figure(5), clf;
FS = 14; % fontsize

% Analyse all datasets, or load dataset
DoAnalyses = 1;
if DoAnalyses == 1
    CPI = []; CPC = [];
    for DS = 1:111
        disp(['Starting dataset #' num2str(DS)])
        % grab fitted data set and calcualte data
        if DS == 111
            [cpi,cpc] = sub_MEP(squeeze(FittingResults(DS,:,:)),1); % also plot
        else
            [cpi,cpc] = sub_MEP(squeeze(FittingResults(DS,:,:)),0);
        end
        
        % Store the data, ignore diagonals (self-impacts)
        m = length(cpi);
        cpi = cpi(:); cpi(1:m+1:end) = [];
        cpc = cpc(:); cpc(1:m+1:end) = [];
        CPI = [CPI; cpi(:)]; CPC = [CPC; cpc(:)];
    end
    % save results
    save DoneAnalysesPerturb
else
    % load results
    load DoneAnalysesPerturb
    sub_MEP(squeeze(FittingResults(108,:,:)),1); % create plots for specific network
end

% plot results for all networks
load SharedParameters CMP 
DefinitionAmbiguous = [0.15 0.17];

Edges = linspace(0,1,25);
for sp = 1:2
    figure(3+sp), hold on
    subplot(2,5,[4,5,9,10]), cla,hold on
    if sp == 1
        H = histc(CPC,Edges); H = H./sum(H); M = ones(size(Edges));
        set(gca,'xtick',[0 0.15 0.4 0.6 0.85 1],'fontsize',FS)
        ylabel('Frequency','fontsize',FS)
        xlabel('Proportion of models','fontsize',FS)
    else
        H = histc(CPI,Edges); H = H./sum(H); M = ones(size(Edges));
        set(gca,'xtick',[0 0.15 0.4 0.6 0.85 1],'fontsize',FS)
        ylabel('Frequency','fontsize',FS)
        xlabel('Proportion of models','fontsize',FS)
    end 
    B = bar(Edges,H.*M,1); set(B,'edgecolor','none','facecolor',CMP(1,:))
    M = ones(size(Edges)); H(Edges<DefinitionAmbiguous(1) | Edges>1-DefinitionAmbiguous(1)) = 0;
    B = bar(Edges,H.*M,1); set(B,'edgecolor','none','facecolor',CMP(3,:))
end
% save figure
% Make_TIFF('Figures/Figure_Predictions.tiff',[0 0 35 12])

function [Consistent_predict_improve,Consistent_predict_change] = sub_MEP(Fits,PT)
warning('off','all')

% delete the worst fits
load SharedParameters Q_threshold
SSD = [Fits{:,2}];
Fits = Fits(1:length(SSD),:);
F = find(SSD < quantile(SSD,Q_threshold));
Fits(F,:) = [];

% How many species are there in this experiment
NumSpp = length(Fits{1,4});

% How much are we going to perturb the species by?
Action_amount = 2;
Action_time = 0.5;

T_obs = linspace(0,19,100);
T_post = linspace(0,Action_time,Action_time*20);

figure(6),clf, hold on
% Solve the ODEs given the fitted parameter values
for Action_species = 1:NumSpp;
    
    for pv = 1:size(Fits,1)
        
        % store parameters
        Thisni = Fits{pv,3};
        Thisr =  Fits{pv,4};
        ThisA =  Fits{pv,5};
        
        % Solve for the observation period
        [t_pv_obs{pv},n_pv_obs{pv}] = ode45(@(t,n)species_DE(t,n,ThisA,Thisr),T_obs,Thisni);
        
        % Solve WITH a press perturbation
        [t_pv_action{pv},n_pv_action{pv}] = ode45(@(t,n)species_DE(t,n,ThisA,Thisr,Action_amount,Action_species),T_post,n_pv_obs{pv}(end,:));
        
        % Solve WITHOUT a press perturbation
        [t_pv_inaction{pv},n_pv_inaction{pv}] = ode45(@(t,n)species_DE(t,n,ThisA,Thisr,0,Action_species),T_post,n_pv_obs{pv}(end,:));
        
        % Does the population increase with the press perturbation?
        Pop_change(pv,:) = n_pv_action{pv}(end,:) > n_pv_action{pv}(1,:);
        
        % Are outcomes better with the press perturbation (i.e., compared to doing nothing)?
        Improved_outcome(pv,:) = n_pv_action{pv}(end,:) > n_pv_inaction{pv}(end,:);
        
        % clear
        clearvars -except *ction* T_* *pv* *_search Tr nr Fits Pop_change Improved_outcome Consistent* NumSpp FileName PT
    end
    
    % store outcomes
    Consistent_predict_improve(Action_species,:) = mean(Improved_outcome);
    Consistent_predict_change(Action_species,:) = mean(Pop_change);
end

if PT == 1
    % Create a colormap where blue is ambiguous, and red is consistent
    load SharedParameters DefinitionAmbiguous CMP
    DefinitionAmbiguous = [0.15 0.17];
    colormap(CMP), set(gcf,'color','w'); FS = 14;
    figure(4),hold on, subplot(2,5,[1 2 6 7]); 
    
    % plot if population increases with perturbation
    pcolor_mike(Consistent_predict_change,1:NumSpp,1:NumSpp,'none',CMP,DefinitionAmbiguous)
    YT = [0.03 0.07 0.14 0.19 0.5]; YT = [YT 1-YT(4:-1:1)];
    YTL = {'decrease';'Definite';'decrease';'Likely';'Ambiguous';'increase';'Likely';'increase';'Definite'};
    box on, axis equal tight off
    for i = 1:NumSpp
        text(i+0.5,0.9,['$s_' num2str(i) '$'],'fontsize',FS,'horizontalalignment','center','Interpreter','latex')
        text(0.9,i+0.5,['$s_' num2str(i) '$'],'fontsize',FS,'horizontalalignment','center','Interpreter','latex')
    end
    text(0.5,3.5,'Perturbed species','fontsize',FS,'rotation',90,'horizontalalignment','center','Interpreter','latex')
    text(3.5,0.5,'Response species','fontsize',FS,'horizontalalignment','center','Interpreter','latex')
    title('Increase or decrease?','fontsize',FS+2,'Interpreter','latex')
    xlim([0.9 NumSpp+1.1])
    ylim([0.9 NumSpp+1.1])
    
    % plot if population is better with perturbation
    figure(5),hold on,subplot(2,5,[1 2 6 7]);
    pcolor_mike(Consistent_predict_improve,1:NumSpp,1:NumSpp,'none',CMP,DefinitionAmbiguous);
    box on, axis equal tight off
    for i = 1:NumSpp
        text(i+0.5,0.9,['$s_' num2str(i) '$'],'fontsize',FS,'horizontalalignment','center','Interpreter','latex')
        text(0.9,i+0.5,['$s_' num2str(i) '$'],'fontsize',FS,'horizontalalignment','center','Interpreter','latex')
    end
    text(0.5,3.5,'Perturbed species','fontsize',FS,'rotation',90,'horizontalalignment','center','Interpreter','latex')
    text(3.5,0.5,'Response species','fontsize',FS,'horizontalalignment','center','Interpreter','latex')
    title('Better or worse?','fontsize',FS+2,'Interpreter','latex')
    xlim([0.9 NumSpp+1.1])
    ylim([0.9 NumSpp+1.1])
end

% save figure
% Make_TIFF(['Experimental fitting results/AmbiguityMatrix_' FileName],[0 0 15 25],'-r100')
return

% Plot the timeseries, before and after the press perturbation
figure(5), clf; CL = get(gca,'colororder'); set(gcf,'color','w'); FS = 14;
for s = 1:4
    subplot(2,2,s), hold on, box on
    if s == Action_species; title('This species was controlled','fontsize',FS,'color','b'); set(gca,'xcolor','b','ycolor','b'); end
    for pv = 1:length(V_fit_search)
        plot(t_pv_obs{pv},n_pv_obs{pv}(:,s),'-','linewidth',1,'color',0.7.*ones(1,3))
        plot([t_pv_obs{pv}; t_pv_obs{pv}(end)+t_pv_action{pv}],[n_pv_obs{pv}(:,s); n_pv_action{pv}(:,s)],'-','linewidth',1,'color',0.7.*ones(1,3))
    end
    plot(Tr,nr(:,s),'.','markersize',18,'color',CL(s,:))
    plot([20 20],[0 4],'r--')
    xlabel('Timesteps elapsed','fontsize',FS); ylabel('Abundance','fontsize',FS);
    ylim([0 1.1.*max(nr(:,s))]); xlim([0 t_pv_obs{pv}(end)+t_pv_action{pv}(end)])
end
% save figure
% Make_TIFF(['Experimental fitting results/AbundanceDynamics_' FileName],[0 0 15 25],'-r100')
