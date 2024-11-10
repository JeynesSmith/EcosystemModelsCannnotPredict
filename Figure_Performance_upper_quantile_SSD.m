% This script creates Supplementary Figure 1. This script asseses a select
% number of experimental datasets, and the ensemble of model fits which
% were found. Using standard fiting methods we identified a series of local
% minima which were believed to be sufficiently good fits according to the
% fitting procedure. For all parameter fits, we display the distribution 
% of their associated Sum of Squared deviations (SSD) (fit statistic). We 
% then highlight the best 10% of fits and display the distribution of 
% their SSDs.

clear all

load 'Experimental fitting results/Results' FittingResults

figure(1), clf, hold on; FS = 18;

Q = 0.1; c = 0;
DS = [111 65 53 59 randi(111,[1,3])]
for c = 1:4
    ds = DS(c);

    SSD = [FittingResults{ds,:,2}];
    qq = quantile(SSD,Q);
    mm = min(SSD);
    rel(ds) = qq./mm - 1;

    edges = linspace(min(SSD),max(SSD),50);

    subplot(2,2,c), hold on
    set(gca,'fontsize',FS-5)

    H = histc(SSD,edges);
    b = bar(edges,H./length(SSD),1);
    set(b,'edgecolor','none')

    H = histc(SSD(SSD<qq),edges);
    b = bar(edges,H./length(SSD),1);
    set(b,'edgecolor','none')

    if c == 1
        L = legend('All local minima','Equivalently good fits');
        set(L,'box','off','fontsize',FS)
    end

    xlabel('SSD of calibrated model','fontsize',FS)
    ylabel('Frequency','fontsize',FS)
end

% Make_TIFF('Figure_S1.tiff',[0 0 28 20])