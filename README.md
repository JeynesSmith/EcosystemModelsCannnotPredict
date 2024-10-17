# EcosystemModelsCannnotPredict
Ecosystem models cannot predict the consequences of conservation decisions

This repository contains the data and computer code (in Matlab) for all mathematical analysis, and for generating all figures, in the publication "Ecosystem models cannot predict the consequences of conservation decisions", submitted to Proceedings of the National Academy of Sciences.

Authors: Larissa Lubiana Botelho, Cailan Jeynes Smith, Sarah Vollert, Michael Bode, December 2023.

The data used in this project was collected by Pennekamp et al. (2018), and can be found in a Matlab compatible format, ‘Analysis_Timeseries.mat’. The experimental data is then fit by running ‘Fit_all_experimental_data.m’

The following scripts were used to generate figures for the above-noted publication. Specifically:

Run SetSharedParameters.m for all figures.
Figure 1 uses Figure_Example_Fits.m
Figure 2 uses Figure_parameters_vs_predictions_systematic.m
Figure 3 uses Figure_predict_responses_to_perturbation_Split.m
Figure 4 uses Figure_MechanisticsAmbiguitySquare.m
S.Figure 1 uses Figure_Performance_upper_quantile_SSD.m
S.Figures 2-5 use Figure_Predictive_accuracy.m
S.Figures 6,7 use Predictive_accuracy_Alternative_Null.m
S.Figures 8-22 use parameter sets fits with alternative fitting procedures and 'ProposedTargets_Unique.mat'.

For any further questions, please contact larissa.lubianabotelho@hdr.qut.edu.au.
