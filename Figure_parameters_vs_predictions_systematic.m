% This script systematically searches through all experimental datasets for 
% parameter fits to the dataset which: a) have highly correlated parameters 
% that produce vastly different predictions, and b) uncorrelated/negatively 
% correlated parameter sets which produce near-identical predictions.

clear all

for reps = 1:10 % loop of 10 iterations of each
for DS = 1:111 % loop of data sets
    Performance = Search_for_system(DS)

    if Performance(1,1) > 0 & Performance(2,1) < 0 % if parameters have been found that meet conditions, then save the figure
        FileName = ['Figures/Comps/Figure_Correlation_ds_' num2str(DS) ...
            '_ps1_' num2str(Performance(1,2)) ...
            '_' num2str(Performance(1,3)) ...
            '_ps2_' num2str(Performance(2,2)) ...
            '_' num2str(Performance(2,3)) ...
            '_rep_' num2str(reps) ...
            '.TIFF']
        Make_TIFF(FileName,[0 0 35 20])
        
    end

end
end



function Performance = Search_for_system(DS)

warning off
% load fits and experimental data
load 'Experimental fitting results/Results' FittingResults NumSets
load 'Experimental data'/Analysis_Timeseries
load SharedParameters

% Call up the data set in question
TS = Analysis_Timeseries{DS,4};

% delete the worst fits
Q_threshold = 0.5;
Fits = squeeze(FittingResults(DS,:,:));
SSD = [Fits{:,2}]; Fits = Fits(1:length(SSD),:);
F = find(SSD > quantile(SSD,Q_threshold));
Fits(F,:) = [];

% define data parameters
NumSpp = length(Fits{1,3}); % number of species
NumSets = size(Fits,1); % number of fits

RHO_thresholds = [0.9 -0.3];

% Set up the figure
figure(3), clf; FS = 16; FSA = 12; MS = 15; LW = 2;
CL = get(gca,'colororder');

CC = [0 0]; CCnum = 0;
% loop over random fits
Performance = zeros(2,3);
for j = randperm(NumSets)%1:(NumSets-1)
    for k = randperm(NumSets)%(j+1):NumSets

        if CCnum == 2
            break
        end

        % Grab two parameter sets and calculate their correlation
        R = [j,k];
        rho_parameters = corrcoef((Fits{R(1),5}(:)),(Fits{R(2),5}(:)));
        rho_parameters = rho_parameters(1,2);

        % only simulate if rho_parameters will meet conditions and
        % parameter sets have not already been identified
        if (rho_parameters < RHO_thresholds(2) && CC(2) == 0) || (rho_parameters > RHO_thresholds(1) && CC(1) == 0)

            % simulate the values for 25 timesteps
            T_obs = linspace(0,25,100);

            % assign network parameters for model 1
            Thisni = Fits{R(1),3};
            Thisr =  Fits{R(1),4};
            ThisA =  Fits{R(1),5};
            % simlate model 1
            [t1,n1] = ode45(@(t,n)species_DE(t,n,ThisA,Thisr),T_obs,Thisni);
            % assign network parameters for model 2
            Thisni = Fits{R(2),3};
            Thisr =  Fits{R(2),4};
            ThisA =  Fits{R(2),5};
            % simlate model 2
            [t2,n2] = ode45(@(t,n)species_DE(t,n,ThisA,Thisr),T_obs,Thisni);

            % scale abundances
            for i = 1:NumSpp
                N = max(TS(i,:));
                n1(:,i) = n1(:,i)./N;
                n2(:,i) = n2(:,i)./N;
            end

            % check that simulated abundances are correct length
            if size(n1,1)==100 && size(n2,1)==100
                % calculate similarity between final abundance of models
                SSD = sum((n1(end,:)-n2(end,:)).^2);
                % calculate similarity between timeseries of models
                fitSSD = sum(sum((n1(1:75,:)-n2(1:75,:)).^2));
                if fitSSD > 0.85; break; end

                % Plot two timeseries with the similar parameters but different predictions
                if SSD > 0.5 && rho_parameters > RHO_thresholds(1) && CC(1) == 0 && length(n1)==length(n2)

                    % plot correlation between parameters
                    if NumSpp == 3
                        subplot(NumSpp,3,2)
                    else
                        subplot(NumSpp,3,[2 1+NumSpp])
                    end

                    plot((Fits{R(1),5}(:)),(Fits{R(2),5}(:)),'.','markersize',MS+8,'color',CL(1,:)), axis equal, axis square
                    hold on, plot([min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))],[min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))],'k--')
                    xlabel('$\alpha_{ij}$, Model 1','fontsize',FS,'Interpreter','latex')
                    ylabel('$\alpha_{ij}$, Model 2','fontsize',FS,'Interpreter','latex')
                    xlim([min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))])
                    ylim([min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))])

                    % plot timeseries
                    for i = 1:NumSpp
                        subplot(NumSpp,3,3*(i-1)+1), hold on, box on
                        plot(t1,n1(:,i),'-' ,'color',CL(1,:),'linewidth',1.5)
                        plot(t1,n2(:,i),'--','color',CL(1,:),'linewidth',1.5)
                        plot([0:18],TS(i,:)./max(TS(i,:)),'o','markersize',5,'color','k','linewidth',1.5)
                        plot([18.5 18.5],[0 2],'k:')
                        ylim([0 2]); xlim([-0.25 max(T_obs)+0.5])
                        ylabel(['$n_{' num2str(i) '}(t)$'],'fontsize',FS,'Interpreter','latex')
                    end
                    xlabel('Timesteps','fontsize',FS,'Interpreter','latex')
                    goodRho=sum(sum((n1(1:75,:)-n2(1:75,:)).^2))
                    CCnum = CCnum+1; CC(1)=1;
                    Performance(1,1) = rho_parameters; % parameter correlation
                    Performance(1,2) = j; % parameter set 1
                    Performance(1,3) = k; % parameter set 2

                end

                % Plot two timeseries with different parameters but the same predictions
                if SSD < 0.05 && rho_parameters < RHO_thresholds(2) && CC(2) == 0 && length(n1)==length(n2)
                    % plot correlation between parameters
                    if NumSpp == 3
                        subplot(NumSpp,3,8)
                    else
                        subplot(NumSpp,3,[3*(NumSpp-1)-1 3*NumSpp-1])
                    end
                    plot((Fits{R(1),5}(:)),(Fits{R(2),5}(:)),'.','markersize',MS+8,'color',CL(2,:)), axis equal, axis square
                    hold on, plot([min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))],[min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))],'k--')
                    xlabel('$\alpha_{ij}$, Model 1','fontsize',FS,'Interpreter','latex')
                    ylabel('$\alpha_{ij}$, Model 2','fontsize',FS,'Interpreter','latex')
                    xlim([min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))])
                    ylim([min(min(([Fits{R(1),5}(:) Fits{R(2),5}(:)]))) max(max(([Fits{R(1),5}(:) Fits{R(2),5}(:)])))])

                    % plot timeseries
                    for i = 1:NumSpp
                        subplot(NumSpp,3,3*i), hold on, box on
                        plot(t1,n1(:,i),'-' ,'color',CL(2,:),'linewidth',1.5)
                        plot(t1,n2(:,i),'--','color',CL(2,:),'linewidth',1.5)
                        plot([0:18],TS(i,:)./max(TS(i,:)),'o','markersize',5,'color','k','linewidth',1.5)
                        plot([18.5 18.5],[0 2],'k:')
                        ylim([0 2]); xlim([-0.25 max(T_obs)+0.5])
                        ylabel(['$n_{' num2str(i) '}(t)$'],'fontsize',FS,'Interpreter','latex')
                    end
                    xlabel('Timesteps','fontsize',FS,'Interpreter','latex')
                    badRho = sum(sum((n1(1:75,:)-n2(1:75,:)).^2))
                    CCnum = CCnum+1; CC(2)=1;

                    Performance(2,1) = rho_parameters; % parameter correlation
                    Performance(2,2) = j; % parameter set 1
                    Performance(2,3) = k; % parameter set 2
                end
            end
        end
    end
end
end










