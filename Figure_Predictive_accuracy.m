% This script is used to analyse the predictive accuracy of models fit to
% all experimental datasets. This information is then displayed in
% Supplementary Figures 2-5. 


%% find fits and get observed + predicted values
compare.data = [];
compare.model = [];
Target = 500; % How many fits should we target?
DS_max = 111; % Number of datasets (experiments)
n_keep = 100; %number of best fit models to keep

errorMargin =  [0.5 0.25 0.1]; %(e.g. 50% error margin) must be ordered least to most constraining
no_errorMargins = size(errorMargin,2);
n_species = zeros(1,DS_max);

% Define indicators for whether all species predictions are within the error
% margins
inside_allspecies = zeros(no_errorMargins,DS_max); %number of models predicting inside the error margin
inside_allspecies_percent = zeros(no_errorMargins,DS_max); %percentage of models predicting inside the error margin
closer_allspecies = zeros(1,DS_max); %indicator on how many times model is closer to predicting true abundance (compared to null)
closer_allspecies_percent = zeros(1,DS_max); %percentage of times the model is closer to predicting true abundance (compared to null)

%Define indicators for whether a single species prediction is within the
%error margins
inside_onespecies = zeros(no_errorMargins,DS_max);
inside_onespecies_percent = zeros(no_errorMargins,DS_max);
closer_onespecies = zeros(1,DS_max);
closer_onespecies_percent = zeros(1,DS_max);

rels_null = [];
mean_rels = [];

%Experiment to plot to view timeseries forecasts
exp_to_plot = 91;

for DS = 1:DS_max
    fprintf('Running experiment %d\n',DS)

    %% Keep only the best fits
    %load all fits
    load('Experimental fitting results/Results','FittingResults')

    %calculate and sort SSDs
    SSDs = [FittingResults{DS,:,2}];
    [~, order] = sort(SSDs);

    %only keep fits with best SSD
    BestFits = squeeze(FittingResults(DS,order(1:n_keep),:));

    %% Get timeseries info for this experiment
    %get full timeseries
    load 'Experimental data'/Analysis_Timeseries.mat
    TS_full = Analysis_Timeseries{DS,4};

    Tr_full = [1:size(TS_full,2)]'; %times
    nr_full = TS_full'; %observations
    n_species(DS) = size(BestFits{1,3},1);

    %% Calc relative error
    %load model predictions and observations
    model_prediction = zeros(n_keep,n_species(DS));
    observation_real = model_prediction;
    for i=1:n_keep
        model_prediction(i,:) = BestFits{i,7}(end,:);
        observation_real(i,:) = nr_full(end,:);
    end

    %calculate relative error
    relative_err = (model_prediction-observation_real)./observation_real; %for each GLV model
    mean_rels = [mean_rels mean(relative_err)]; %average relative error across each GLV model
    rels_null = [rels_null (nr_full(end-1,:)-observation_real(1,:))./observation_real(1,:)]; %for null model

    %absolute error
    abs_err = (model_prediction-observation_real);

    %% print the acceptable error

    %if this this the experiment that we with to plot
    if exp_to_plot == DS
        %plot
        figure;
        CL = [0.4627    0.7098    0.7725;0.1176    0.5059    0.6902; 0.0824    0.2980    0.4745];
        MXX = max(nr_full); FS = 11;
        tiles = tiledlayout(n_species(DS),1);
    end

    %for each error margin specified
    for er=1:no_errorMargins
        %define the range of predictions within the error margin
        upper_bound = observation_real*(1+errorMargin(er));
        lower_bound = max(observation_real*(1-errorMargin(er)),0);

        %PREDICTIONS OF ALL SPECIES
        %for each model
        for i=1:n_keep
            %check if within bounds
            if [lower_bound(i,:)<model_prediction(i,:) model_prediction(i,:)<upper_bound(i,:)]
                %count model as within bounds
                inside_allspecies(er, DS) = inside_allspecies(er, DS) +1;
            end

            %check if it beats the null model
            if abs(model_prediction(i,:)-observation_real(i,:)) < abs(nr_full(end-1,:)-observation_real(i,:)) %if LV is closer
                %count model as beating the null model
                closer_allspecies(DS) = closer_allspecies(DS)+1;
            end
        end
        %calculate the percentages of models inside error margins and beating null model
        inside_allspecies_percent(er, DS) = inside_allspecies(er, DS)/n_keep;
        closer_allspecies_percent(DS) = closer_allspecies(DS)/n_keep;

        %PREDICTIONS OF INDIVIDUAL SPECIES
        for i=1:n_keep         %for each model

            for s=1:n_species(DS)       %for each species

                %check if prediction is within bounds
                if [lower_bound(i,s)<model_prediction(i,s) model_prediction(i,s)<upper_bound(i,s)]
                    %count model as within bounds
                    inside_onespecies(er, DS) = inside_onespecies(er, DS) +1;
                end

                %check if it beats the null model
                if abs(model_prediction(i,s)-observation_real(i,s)) < abs(nr_full(end-1,s)-observation_real(i,s)) %if LV is closer
                    %count model as beating the null model
                    closer_onespecies(DS) = closer_onespecies(DS)+1;
                end
            end
        end
        %calculate the percentages of models inside error margins and beating null model
        inside_onespecies_percent(er, DS) = inside_onespecies(er, DS)/(n_keep*n_species(DS));
        closer_onespecies_percent(DS) = closer_onespecies(DS)/(n_keep*n_species(DS));


        if exp_to_plot == DS
            % PLOTTING DATA & PREDICTIONS

            %plot outside error margins first
            if er ==1
                for s = 1:n_species(DS)
                    nexttile(s)
                    hold on
                    for i=1:n_keep
                        if (lower_bound(i,s)>model_prediction(i,s) || model_prediction(i,s)>upper_bound(i,s)) %outside acceptable bounds
                            hold on
                            %model fit to data (past)
                            plot(BestFits{i,6}(1:end-5),BestFits{i,7}(1:end-5,s),'-','linewidth',1,'color',[.7 .7 .7]) %plot in grey
                            %model prediction (future)
                            plot(BestFits{i,6}(end-5:end),BestFits{i,7}(end-5:end,s),'-','linewidth',1,'color',[.7 .7 .7]) %plot in grey
                        end
                    end
                end
            end

            %plot inside error margins
            for s = 1:n_species(DS)
                nexttile(s)
                hold on
                %plot all preditions
                for i=1:n_keep
                    if lower_bound(i,s)<model_prediction(i,s) && model_prediction(i,s)<upper_bound(i,s) %within acceptable bounds
                        hold on
                        %model fit to data (past)
                        plot(BestFits{i,6}(1:end-5),BestFits{i,7}(1:end-5,s),'-','linewidth',1,'color',CL(er,:)) %plot in colour       %past
                        %model prediction (future)
                        plot(BestFits{i,6}(end-5:end),BestFits{i,7}(end-5:end,s),'-','linewidth',1,'color',CL(er,:)) %plot in colour       %future
                    end
                end
            end 
        end % plot

    end %error margin


    %plot data
    if exp_to_plot == DS

        for s = 1:n_species(DS)
            nexttile(s)
            hold on
            plot(Tr_full,nr_full(:,s),'.','markersize',18,'color','k')
            xline(18,'--','LineWidth',1.5)
            ylabel(['$n_{' num2str(s) '}(t)$'],'fontsize',FS,'Interpreter','latex')
            xlim([0 19])
            ylim([0 MXX(s)*1.3])
        end
        xlabel('Timesteps','fontsize',FS,'Interpreter','latex')
    end

end


%% Plot histogram for being within error bounds 
CL = [0.4627    0.7098    0.7725;0.1176    0.5059    0.6902; 0.0824    0.2980    0.4745];

%one species at a time
figure; hold on;
tiles = tiledlayout(1,length(errorMargin));
for i=1:no_errorMargins %for each error margin
    nexttile(i)
    histogram(inside_onespecies_percent(i,:),"BinWidth",0.1,"FaceColor",CL(i,:))
    lable_entry(i) = {strcat("$\pm",num2str(round(errorMargin(i)*100))+"$\% error margin")};
    xlim([0 1])
    xlabel(lable_entry(i),'Interpreter',"latex")
end
ylabel(tiles,'Number of experiments','Interpreter',"latex")
xlabel(tiles,'Individual species predictions within error margin of true observation','Interpreter',"latex")


%all species simultaneously
figure; hold on
for i=1:no_errorMargins
    histogram(inside_allspecies_percent(i,:),"BinWidth",0.02,"FaceColor",CL(i,:))
    legend_entry(i) = {strcat("$\pm",num2str(round(errorMargin(i)*100))+"$\% error margin")};
end
xlim([0 1])
xlabel('\% models with all species simultaneously predicted within error margin of true observation','Interpreter',"latex")
ylabel('Number of experiments','Interpreter',"latex")
legend(legend_entry,'Interpreter',"latex")

%% Plot histogram for outperforming null model (can remove)
%outperforming the null model
figure
tiledlayout(1,2)

%one species at a time
nexttile
histogram(closer_onespecies_percent,"BinWidth",0.1)
xlim([0 1])
label = strcat("\% Lokta Volterra model predictions which outperformed the null model");
xlabel(label,'Interpreter',"latex")
ylabel('Number of experimental ecosystems','Interpreter',"latex")

%all species simultaneously
nexttile
histogram(closer_allspecies_percent,"BinWidth",0.1)
xlim([0 1])
label = strcat("\% Lotka Volterra models which outperformed the null model for all species predictions simultaneously");
xlabel(label,'Interpreter',"latex")
ylabel('Number of experimental ecosystems','Interpreter',"latex")


%% Violin plot of prediction performance null model vs GLV
figure; hold on;
cats = categorical({'Lotka-Volterra','Null model'});
violin = violinplot([log(abs(mean_rels)); log(abs(rels_null))]',cats, ...
    'ShowData', false, 'ShowMean',false,'ShowBox',false,'ShowMedian',true,'ShowWhiskers',false);
ylabel('Logarithmic relative error')


%% Show an example of GLV prediction vs null model
DS = exp_to_plot;

%set up figure
figure, clf; CL = get(gca,'colororder'); MXX = max(nr_full)*2; FS = 11;
MXX = [5, 3, 3, 7.2];
tiles = tiledlayout(n_species(DS),1);

%for each species
for s = 1:n_species(DS)
    nexttile
    hold on
    
    %plot up to t=18
    for j=1:5 %for 5 random models

        i=randsample(n_keep,1);
        plot(BestFits{i,6}(1:ceil(100/19*18)),BestFits{i,7}((1:ceil(100/19*18)),s),'-','linewidth',1,'color',CL(1,:)) %plot in colour

        %plot GLV predition
        plot(BestFits{i,6}(ceil(100/19*18):end),BestFits{i,7}((ceil(100/19*18):end),s),':','linewidth',1,'color',CL(1,:)) %plot in colour
        plot(BestFits{i,6}(end), BestFits{i,7}(end,s),'.','markersize',18,'color',CL(1,:))
    end

    %plot null prediction
    plot([18 19],[nr_full(end-1,s) nr_full(end-1,s)],':','linewidth',1,'color',CL(2,:))
    plot(19, nr_full(end-1,s),'.','markersize',18,'color',CL(2,:))

    %plot datapoints
    plot(Tr_full,nr_full(:,s),'.','markersize',18,'color','k')

    %plot endtime
    xline(18)

    ylim([0 MXX(s)])
    ylabel(['$n_{' num2str(s) '}(t)$'],'fontsize',FS,'Interpreter','latex')
    xlim([0 19])
end
xlabel('Timesteps','fontsize',FS,'Interpreter','latex')




