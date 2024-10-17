% This script is responsible for generating Supplementary Figures 6,7 on
% alternate null models. Null models are compared to generated data fits.


%% find fits and get observed + predicted values
SM = 0; % Should we use the smoothed datasets?
compare.data = [];
compare.model = [];
Target = 500; % How many fits should we target?
DS_max = 111;
n_keep = 100; %number of best fit models to keep

rels_nulls = [];
rels_GLV = [];
abs_nulls = [];
abs_GLV = [];

%Experiment to plot to view timeseries forecasts
exp_to_plot = 71;

for DS = 1:DS_max
    fprintf('Running experiment %d\n',DS)

    sub_Multiple_Experimental_Fits(DS,Target);

    %% Keep only the best fits

    %load all fits
    load 'Experimental fitting results'\Results FittingResults %%%%%%%%%%%%%%%%%%%%%%%%%%

    %calculate and sort SSDs
    SSDs = [FittingResults{DS,:,2}];
    [~, order] = sort(SSDs);

    %only keep the best fit
    BestFits = squeeze(FittingResults(DS,order(1),:));

    %% Get timeseries info for this experiment
    %get full timeseries
    load 'Experimental data'\Analysis_Timeseries
    if SM == 0      TS_full = Analysis_Timeseries{DS,4};
    elseif SM == 1  TS_full = Analysis_Timeseries{DS,5};
    end

    Tr_full = [1:size(TS_full,2)]'; %times
    nr_full = TS_full'; %observations
    n_species(DS) = size(BestFits{3},1);

    %% Get predictions for each model
    clearvars observation_real null_predictions GLV_prediction;
    
    %GLV
    GLV_predictions = BestFits{7}(end,:);

    %true observation
    observation_real = nr_full(end,:);

    %null models
    null_predictions = zeros(4,n_species(DS));

    %0th order null model
    x0 = ones(1,n_species(DS))*0.1;
    [x0_fit, ~] = fmincon(@(x) Fit_Null_0th(x, Tr_full(1:end-1),nr_full(1:end-1,:)),x0,[],[],[],[]);
    null_predictions(1,:) = x0_fit;

    %1st order null model
    x1 = ones(1,n_species(DS)*2)*0.1;
    [x1_fit, ~] = fmincon(@(x) Fit_Null_1st(x, Tr_full(1:end-1),nr_full(1:end-1,:)),x1,[],[],[],[]);
    for s=1:n_species(DS); null_predictions(2,s) = x1_fit(2*s-1)*19+x1_fit(2*s); end

    %2nd order null model
    x2 = ones(1,size(nr_full,2)*3)*0.1;
    [x2_fit, ~] = fmincon(@(x) Fit_Null_2nd(x, Tr_full(1:end-1),nr_full(1:end-1,:)),x2,[],[],[],[]);
    for s=1:n_species(DS); null_predictions(3,s) = x2_fit(3*s-2)*19^2+x2_fit(3*s-1)*19+x2_fit(3*s); end

    %exponential model
    x_exp = ones(1,size(nr_full,2)*2)*0.1;
    [x_exp_fit, ~] = fmincon(@(x) Fit_Null_Exp(x, Tr_full(1:end-1),nr_full(1:end-1,:)),x_exp,[],[],[],[]);
    for s=1:n_species(DS); null_predictions(4,s) = x_exp_fit(s*2-1)^19*x_exp_fit(s*2); end
   
    %% Calc relative error

    %calculate relative error
    rels_GLV = [rels_GLV (GLV_predictions-observation_real)./observation_real]; % relative error for best GLV model                                                   %%%%%%%%%%
    rels_nulls = [rels_nulls (null_predictions-observation_real)./observation_real]; %for null models                                                   %%%%%%%%%%

    abs_GLV = [abs_GLV abs(GLV_predictions-observation_real)]; %absolute error for GLV                                                %%%%%%%%%%
    abs_nulls = [abs_nulls abs(null_predictions-observation_real)]; %for null models                                                   %%%%%%%%%%

    %% Plot a comparison
    if rand<0.1 %plot 10% of figures randomly
        %set up figure
        figure, clf; CL = get(gca,'colororder'); MXX = max(nr_full)*1.1; FS = 11;
        tiles = tiledlayout(n_species(DS),1);

        %for each species
        for s = 1:n_species(DS)
            nexttile
            hold on

            %plot GLV model fit and prediction
            plot(BestFits{6}(1:ceil(100/19*18)),BestFits{7}(1:ceil(100/19*18),s),'-','linewidth',1,'color',CL(1,:))
            plot(BestFits{6}(ceil(100/19*18):end),BestFits{7}(ceil(100/19*18):end,s),':','linewidth',1,'color',CL(1,:)) %plot in colour
            plot(BestFits{6}(end), BestFits{7}(end,s),'.','markersize',18,'color',CL(1,:))

            %plot null predictions
            %0th order
            null_0 = @(s) x0_fit(s);
            plot([0 18],[null_0(s) null_0(s)],'-','linewidth',1,'color',CL(2,:))
            plot([18 19],[null_0(s) null_0(s)],':','linewidth',1,'color',CL(2,:))
            plot(19, null_0(s),'.','markersize',18,'color',CL(2,:))
            %1st order
            null_1 = @(s,t) x1_fit(s*2-1)*t+x1_fit(s*2);
            plot([0 18],[null_1(s,0) null_1(s,18)],'-','linewidth',1,'color',CL(3,:))
            plot([18 19],[null_1(s,18) null_1(s,19)],':','linewidth',1,'color',CL(3,:))
            plot(19, null_1(s,19),'.','markersize',18,'color',CL(3,:))
            %2nd order
            null_2 = @(s,t) x2_fit(s*3-2)*t.^2+x2_fit(s*3-1)*t + x2_fit(s*3);
            plot(0:18,null_2(s,0:18),'-','linewidth',1,'color',CL(4,:))
            plot([18 19],[null_2(s,18) null_2(s,19)],':','linewidth',1,'color',CL(4,:))
            plot(19, null_2(s,19),'.','markersize',18,'color',CL(4,:))
            
            %exponential
            null_5 = @(s,t) x_exp_fit(s*2-1).^t*x_exp_fit(s*2);
            plot(0:18,null_5(s,0:18),'-','linewidth',1,'color',CL(5,:))
            plot([18 19],[null_5(s,18) null_5(s,19)],':','linewidth',1,'color',CL(7,:))
            plot(19, null_5(s,19),'.','markersize',18,'color',CL(5,:))

            %plot datapoints
            plot(Tr_full,nr_full(:,s),'.','markersize',18,'color','k')

            %plot endtime
            xline(18)

            ylabel(['$n_{' num2str(s) '}(t)$'],'fontsize',FS,'Interpreter','latex')
            xlim([0 19])
        end
        xlabel('Timesteps','fontsize',FS,'Interpreter','latex')
    end
    pause(1)
end

%% Violin plot of prediction performance null model vs GLV
figure; hold on;
cats = categorical({'Lotka-Volterra','0th order','1st order','2nd order','Exponential'});
violin = violinplot([log(abs(rels_GLV)); log(abs(rels_nulls(1:4,:)))]',cats, ...
    'ShowData', false, 'ShowMean',false,'ShowBox',false,'ShowMedian',true,'ShowWhiskers',false);
ylabel('Logarithmic relative error','Interpreter','latex')
xlabel('Model','Interpreter','latex')
set(gca,'XTickLabel',cats,'TickLabelInterpreter','latex','XTickLabelRotation',90)


%% Pairwise comparison
figure;yline(0.5);%yline(0.3);yline(0.7)
percent_GLV_beats_null = sum(((abs([rels_nulls(1:4,:)]) - abs(rels_GLV))>0)')'/size(rels_nulls,2)
names = categorical({'0th order','1st order','2nd order','Exponential'});
b = bar(names,percent_GLV_beats_null);
b.FaceColor = 'flat';
b.CData(1:4,:) = CL(2:5,:);
set(gca,'XTickLabel',names,'TickLabelInterpreter','latex','XTickLabelRotation',90)
ylim([0 1])
ylabel('\% of predictions where Lotka Volterra outperforms the alternate model','Interpreter','latex')
xlabel('Model','Interpreter','latex')


% % (abs(rels_nulls) - abs(rels_GLV))>0 %if positive, then the GLV performs better
% figure
% percent_GLV_beats_null = sum(((abs(rels_nulls) - abs(rels_GLV))>0)')'/size(rels_nulls,2)
% names = categorical({'0th order','1st order','2nd order','3rd order','4th order','Exponential'});
% b = bar(names,percent_GLV_beats_null);
% b.FaceColor = 'flat';
% b.CData(1:6,:) = CL(2:7,:);
% ylim([0 1])
% yline(0.5);yline(0.3);yline(0.7)
% ylabel('% of times that Lotka Volterra beat the null model')
% xlabel('Null model')

%Split into size of systems
percent_GLV_beats_null_3 = sum(((abs(rels_nulls(:,1:51*3)) - abs(rels_GLV(1:51*3)))>0)')'/size(rels_GLV(1:51*3),2)
percent_GLV_beats_null_4 = sum(((abs(rels_nulls(:,51*3+1:51*3+49*4)) - abs(rels_GLV(51*3+1:51*3+49*4)))>0)')'/size(rels_GLV(51*3+1:51*3+49*4),2)
percent_GLV_beats_null_5 = sum(((abs(rels_nulls(:,51*3+49*4+1:end)) - abs(rels_GLV(51*3+49*4+1:end)))>0)')'/size(rels_GLV(51*3+49*4+1:end),2)

figure
names = categorical({'0th order','1st order','2nd order','Exponential'});
bar(names,[percent_GLV_beats_null_3 percent_GLV_beats_null_4 percent_GLV_beats_null_5])
ylim([0 1])
yline(0.5)
legend('3 species', '4 species', '5 species','Lotka-Volterra win/loss line')

function [SSD] = Fit_Null_0th(x, Tr, nr)
SSDs = zeros(1,size(nr,2));
%for each species
for s=1:size(nr,2)
    data = nr(:,s);
    const = x(s);
    null = const;
    SSDs(s) = sum((data-null).^2);
end
SSD = sum(SSDs);
end

function [SSD] = Fit_Null_1st(x, Tr, nr)
SSDs = zeros(1,size(nr,2));
%for each species
for s=1:size(nr,2)
    data = nr(:,s);
    slope = x(2*s-1);
    int = x(2*s);
    null = @(x) slope*x+int;
    SSDs(s) = sum((data-null(Tr)).^2);
end
SSD = sum(SSDs);
end

function [SSD] = Fit_Null_2nd(x, Tr, nr)
SSDs = zeros(1,size(nr,2));
%for each species
for s=1:size(nr,2)
    data = nr(:,s);
    a = x(3*s-2);
    b = x(3*s-1);
    c = x(3*s);
    null = @(x) a*x.^2+b*x+c;
    SSDs(s) = sum((data-null(Tr)).^2);
end
SSD = sum(SSDs);
end

function [SSD] = Fit_Null_Exp(x, Tr, nr)
SSDs = zeros(1,size(nr,2));
%for each species
for s=1:size(nr,2)
    data = nr(:,s);
    a = x(s*2-1);
    b = x(s*2);
    null = @(x) a.^x*b;
    SSDs(s) = sum((data-null(Tr)).^2);
end
SSD = sum(SSDs);
end