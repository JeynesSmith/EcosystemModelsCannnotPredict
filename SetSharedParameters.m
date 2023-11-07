clear all
% This file assigns shared parameters which are used across all 
% files. This dataset includes 1) the upper quantile for SSDs kept/used in
% analysis, 2) The colours used for plotting ambiguity, 3) and the
% thresholds for ambiguity. These are saved as SharedParameters.mat

% What quantile of SSD fit results will we keep?
Q_threshold = 0.25;

% What are the colours of ambiguous and unambiguous
CMP =   [[0.7 1 0.7]; ... % Green = unambiguous
        [0.7 0.7 1]; ... % Blue = halfway
        [0.5 0.5 0.5]]; % Grey = ambiguous

DefinitionAmbiguous = [0.15 0.3]; % the thresholds of ambiguity (15% = unambiguous, 30% = halfway)

% save parameters 
save SharedParameters