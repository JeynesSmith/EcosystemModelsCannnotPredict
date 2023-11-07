% This script creates a plot of how consistent the sign of an interaction
% was amongst fitted models. For every interaction in every experimental
% dataset, we plot a single point which identifies how consistent the sign 
% of the interaction was between species i and species j, and the inverse. 
% A result of 1 indicates that species i always positively affected species
% j, whereas a result of 0, indicates that i always had a negative effect
% on j.

clear all

DoAnalyses = 0;
if DoAnalyses == 1
    % Load all the datasets
    load 'Experimental fitting results'/Results FittingResults NumSets
    load SharedParameters

    % Check interaction sign consistency
    InteractionTypeConsistency = [];
    
    % loop over experimental datasets
    for DS = 1:111
                
        % delete the worst fits
        Fits = squeeze(FittingResults(DS,:,:));
        SSD = [Fits{:,2}];
        Fits = Fits(1:length(SSD),:);
        F = find(SSD < quantile(SSD,Q_threshold));
        Fits(F,:) = [];
        
        % determine the number of species and number of fits
        NumSpp = length(Fits{1,3});
        NumSets = size(Fits,1);
        
        % store interaction matrix
        for NS = 1:NumSets
            A(:,:,NS) = Fits{NS,5};
        end
                
        % Go through all the species interactions one-by-one
        for s1 = 1:NumSpp
            for s2 = s1+1:NumSpp
                % calculcate the percentage of times that aij and aji were
                % in the positive region. 
                InteractionTypeConsistency  = [InteractionTypeConsistency; sum(A(s1,s2,:) > 0)./NumSets , sum(A(s2,s1,:) > 0)./NumSets];
            end
        end
        clearvars -except FittingResults ImpactConsistency ...
            ConsistentPredator InteractionTypeConsistency ...
            InteractionTypeConsistency_null DominantCompetitorConsistency ...
            MeanStd Q_threshold CMP DefinitionAmbiguous
    end
    % save analysis
%     save DoneAnalysesSquare
end

% load analysis
clear all
load DoneAnalysesSquare *onsist* MeanStd Q_threshold CMP DefinitionAmbiguous

% Plot the interaction types as a Ternary plot
Amb_level_1 = DefinitionAmbiguous(1);

figure(1), clf; ex = 0.01; set(gcf,'color','w'); box on
subplot('position',[0.1 0.1 0.8 0.8]); hold on,  
MS = 8; FS = 14;

% Plot the main green square
pp = patch([0 0 1 1],[0 1 1 0],CMP(1,:));
set(pp,'edgecolor','none')

greyVector = [Amb_level_1 0; ...
              Amb_level_1 Amb_level_1; ...
              0 Amb_level_1; ...
              0 1-Amb_level_1 ; ...
              Amb_level_1 1-Amb_level_1 ; ...
              Amb_level_1 1; ...
              1-Amb_level_1 1; ...
              1-Amb_level_1 1-Amb_level_1 ; ...
              1 1-Amb_level_1 ; ...
              1 Amb_level_1 ; ...
              1-Amb_level_1 Amb_level_1 ; ...
              1-Amb_level_1 0];

% Plot the central grey square
pp = patch(greyVector(:,1),greyVector(:,2),0.9.*ones(1,3)); set(pp,'edgecolor','none')

% halfway lines
axis square
plot([0 1],[0.5 0.5],'k--')
plot([0.5 0.5 nan 0.5 0.5],[0 0.82 0.82 0.89 1],'k--')

pp = gca; pp.FontSize = FS;
pp.XTick = [0,1];
pp.YTick = [0,0.15,0.5,0.85,1];
xlabel('Proportion $a_{ij}>0$','Interpreter','latex'),ylabel('Proportion $a_{ji}>0$','Interpreter','latex')

% Label
TX = [0.05   -0.065
    0.06 1.065
    0.94 -0.065
    0.94 1.065
    0.5 0.86];

FS = 16;
text(TX(1,1),TX(1,2),'Competition','fontsize',FS,'color',[0 0.35 0],'verticalalignment','bottom','Interpreter','latex','fontsize',FS)
text(TX(2,1),TX(2,2),'$j$ predates $i$','fontsize',FS,'color',[0 0.35 0],'verticalalignment','top','Interpreter','latex','fontsize',FS)
text(TX(3,1),TX(3,2),'$i$ predates $j$','fontsize',FS,'color',[0 0.35 0],'horizontalalignment','right','verticalalignment','bottom','Interpreter','latex','fontsize',FS)
text(TX(4,1),TX(4,2),'Mutualism','fontsize',FS,'color',[0 0.35 0],'horizontalalignment','right','verticalalignment','top','Interpreter','latex','fontsize',FS)
text(TX(5,1),TX(5,2),'Ambiguous','fontsize',FS,'color','k','horizontalalignment','center','Interpreter','latex','fontsize',FS)

% scatter plot the interaction percentages
plot(InteractionTypeConsistency(:,1),InteractionTypeConsistency(:,2),'k.','markersize',MS)
box on
set(gca,'layer','top')
axis square

% save figure
Make_TIFF('Figures/Figure_mechanistic_ambiguity_square_Final.TIFF',[0 0 20 20])


