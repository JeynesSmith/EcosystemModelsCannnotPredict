# EcosystemModelsCannnotPredict
Ecosystem models cannot predict the consequences of conservation decisions

This repository contains the data and computer code (in Matlab 2024a, Windows 11) for all mathematical analysis, and for generating all figures, in the publication "Calibrated ecosystem models cannot predict the consequences of conservation management decisions", submitted to Ecology Letters.

Authors: Larissa Lubiana Botelho, Cailan Jeynes-Smith, Sarah Vollert, Michael Bode, October 2024.

The data used in this project was collected by Pennekamp et al. (2018), excluding experiments in which species have gone extinct. Data is stored in a Matlab compatible format in ‘Analysis_Timeseries.mat’. ‘Analysis_Timeseries.mat’ is a cell format where rows are different experimental datasets and the columns contain: the number of species in the experiment, the raw data from Pennekamp et al., smoothed data, data normalised by the mean for each species, and data normalised by the mean and smoothed. In this study, we have only used the normalised data in the fourth column. 

To prepare their directory, a user must create three subfolders within the working directory: 'Experimental data' (containing Analysis_Timeseries.mat, and 'ProposedTargets_Unique.mat'), 'Experimental fitting results' (where fits will be saved as 'Results.mat'), and 'Figures' (where figures will be saved).

The experimental data is fit by running ‘Fit_all_experimental_data.m’. The saved fits to each experiment can be found in ''Experimental fitting results'\Results.mat' which contains the 3-dimensional cell, 'FittingResults', and the number of fits to each set of experimental data, 'Target'. The first dimension of 'FittingResults' represents the experimental data which the code is fit to, the second dimension represents the individual fits to the dataset, and the third dimension contains: the parameters of the fit, the SSD measurement, the experiment initial conditions, the fit growth rate, the interaction matrix, the time vector for a simulated fit, and the abudnance vector for a simulated fit. We recommend creating at least 500 fits for each experimental dataset. We have provided a small sample of 50 fits for each dataset.

Supplementary fitting procedures are created using 'Supp_Fit_Data_AIC.m', 'Supp_Fit_Data_Bayesian.m', 'Supp_Fit_Data_Lasso.m', and 'Supp_Fit_Data_ThirdOrder_Lasso.m'. These will create the same data structures as 'Fit_all_experimental_data.m', except for 'Supp_Fit_Data_ThirdOrder_Lasso.m' which will have an additional dimension for storing the higher-order term parameters.

This codebase can also be adapted to run in 'GNU Octave Version 8.4.0', an opensource alternative to Matlab. We have provided the adapted code for all primary analysis and main text figures in 'Octave Files'. We recommend using Matlab where possible as this work is not optimised for Octave and therefore requires a significant increase in runtime.

The following scripts were used to generate figures for the above-noted publication. Specifically:

Run 'SetSharedParameters.m' before attempting to generate figures.
Figure 1 uses Figure_Example_Fits.m
Figure 2 uses Figure_parameters_vs_predictions_systematic.m
Figure 3 uses Figure_predict_responses_to_perturbation_Split.m
Figure 4 uses Figure_MechanisticsAmbiguitySquare.m
S.Figure 1 uses Figure_Performance_upper_quantile_SSD.m
S.Figures 2-5 use Figure_Predictive_accuracy.m
S.Figures 6,7 use Predictive_accuracy_Alternative_Null.m
S.Figures 8-22 use parameter sets fits with alternative fitting procedures and 'ProposedTargets_Unique.mat'.

For any further questions, please contact larissa.lubianabotelho@hdr.qut.edu.au.
