% LIP single-trial decision variable analyses


%% General setting
% Always need to run this section

global mLoc
mLoc = 'My\code\location\m-files';
cd(mLoc)
dat_path = 'My\data\location\single trial_data';
global saveLoc 
saveLoc = fullfile(dat_path, 'All Sessions'); 
if exist(saveLoc) == 0
    mkdir(saveLoc)
end
global dNames 
dNames = {'TinC', 'TinI', 'Ramp', 'PC1', 'WhatD', 'WhenD', 'MinC', 'MinI'};

general_settings

%% Preprocessing

preproccess

%% Generate weight vecotrs for identified neuron classes and global parameter file par_sess

comb_par

%% Min indentification

min_identification


%% Compute Ramp direction
% Is working smoothly, but still needs to be commented
disp('Computing ramp directions')
calc_ramp_dim

%% Compute PC1 - Requires mmx
disp('Computing PCA')
calc_pca
calc_pca_Min
% Evaluate PCA results - updated and recomputed!
compute_stats_PCA

%% Compute When-decoder

calc_when_decoders

% calc_when_decoders_choice_predict 

%% Compute What-decoder
disp('Computing decoder directions')

% Caluclate time-dependent What-decoders 
calc_whatDecoder_time_dependent % Last script run during test

% Calculate weight vectors for What-decoder dimension (450ms after motion onset)
calc_whatDecoder_fixed_time

% Calculate decoder performance of across time points (former calc_Fig4b)
calc_whatDecoder_across_tpts         


%% Generate population signals for all sessions
%% Apply weight vectors to population data
%% Combine all signals into one signal structure

sigSL = []; sigRL = []; 
calc_sig_struct

sigSLv = []; sigRLv = [];
calc_sig_struct_view


%% ***********************************************************************
%% FIGURE 1
%% Perceptual decisions are explained by the accumulation of noisy evidence to a stopping bound.
%% ***********************************************************************
% Behavior


%% ***********************************************************************
%% FIGURE 2
%% Population responses from LIP approximate drift-diffusion.
%% ***********************************************************************
% Single decisions and trial averages for Ramp, PC1 and Tin-con - DONE

load_Fig2 

plot_Figure2 


%% ***********************************************************************
%% FIGURE 3
%% Variance and autocorrelation of the single trial signals.
%% ***********************************************************************
% Var and Cor


%% ***********************************************************************
%% FIGURE 4
%% The population signal predictive of choice and RT is approximately one-dimensional.
%% ***********************************************************************

plot_Figure4

calc_Fig4_stats

% includes calc_Fig4de which contains the calculation of stats on CS and
% within-trial correlation (reported in paper)



%% ***********************************************************************
%% FIGURE 5
%% The drift-diffusion signal \rev{approximates the} decision variable.
%% ***********************************************************************
% a) Leverage of activity on choice: Ramp, PC1, Tin-con
% b) Correlation of activity with RT: Ramp, PC1, Tin-con

% Compute mediation
analysis_mediation.run_mediation_all2all

% Plot figure
plot_Figure5

clac_Figure5_stats

%% ***********************************************************************
%% FIGURE 6
%% The representation of momentary evidence in area LIP.
%% ***********************************************************************

plot_Figure6


%% ***********************************************************************
%% SUPPLEMENTARY FIGURES
%% ***********************************************************************


%% ***********************************************************************
%% FIGURE S1
%% Effect of motion pulses on behavior.
%% ***********************************************************************

% Reproduction from Stine et al., 2023. doi: 10.1016/j.neuron.2023.05.028

%% ***********************************************************************
%% FIGURE S2
%% Derivation of a ramp coding-direction in neuronal state space.
%% ***********************************************************************

% plotted in calc_ramp_dim

%% ***********************************************************************
%% FIGURE S3
%% Trial-averaged activity grouped by RT quantile
%% ***********************************************************************
%% Plot motion discrimination task split by RT percentile

figure('Position', [200 10 600 600]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')
plot_FigS3 


%% ***********************************************************************
%% FIGURE S4
%% Trial-averaged activity after subtracting the urgency component
%% ***********************************************************************

plot_FigS4

%% ***********************************************************************
%% FIGURE S5
%% Bounds induce sublinear increase in variance of diffusion paths
%% ***********************************************************************
%% Mike's code

%% ***********************************************************************
%% FIGURE S6
%% Derivation of a When-direction in neuronal state space
%% ***********************************************************************

% plotted in 


%% ***********************************************************************
%% FIGURE S7
%% Single-trials and trial averaged signals furnished by the When- and What-decoders.
%% ***********************************************************************

plot_resp_aligned = 1;
plot_FigS6


%% ***********************************************************************
%% FIGURE S8
%% Variability in cosine similarity and within-trial correlation across sessions.
%% ***********************************************************************

% Cosine similarity & Within-trial correlation for individual sessions
% (Control analyses for Figure 4)
useD = [3 4 1 5 6]; 
dimLabelsCS = {'ramp', 'PC1', 'Tin^{con}_{in}', 'What-decoder', 'When-decoder'};

% 1. Cosine similarity between random projections
calc_weightCS_random
% load(fullfile(saveLoc, 'cosine_sim_randperm'), 'cosine_sim_randperm') % cosine_sim_randperm(sess, d1, d2, :)

% 2. Correlation between single trials in random directions in state space
calc_wiTrialCorr_random

plot_suppFig_CS_trialCorr
plot_suppFig_CS_trialCorr_az


stats_CS_trialCorr_rand


%% ***********************************************************************
%% FIGURE S9
%% The Tin neurons are not discoverable by their weight assignments
%% ***********************************************************************

plot_FigS7


%% ***********************************************************************
%% FIGURE S10
%% Leverage of signals on behavior is non-trivial and as predicted by simulations.
%% ***********************************************************************
%% Control analyses for Figure 5

% 1. Mediation analysis on random projections in the data: op row
mediation_random_projections

% 2. Mediation on simulated diffusion traces: bottom row
% Correlation between Ti nneurons in the data
xneuron_cov
% Mediation on simluated data (correlated neurons)
mediation_simulated

% 3. Mediation with shuffled trial identity (response to reviewer)
mediation_shuffle

% Caclulate mediation statistics
calc_mediation_stats

% Mediation control 
% controlling magnitude of mediation itself by shuffling identity of
% mediating time point
mediation_zeta_control

% Mediation control: individual coherences

mediation_singleCoh


%% ***********************************************************************
%% FIGURE S11
%% Cross-mediation of Single-trial correlations with behavior.
%% ***********************************************************************
% Still need to update all inputs!

calc_FigS8

plot_FigS8


%% ***********************************************************************
%% FIGURE S12
%% Comparison of linear and non-linear choice decoders.
%% ***********************************************************************
% Code from Gabe


%% ***********************************************************************
%% FIGURE S13
%% ***********************************************************************
%% Decoding choice from subsets of neurons

plot_FigS9


%% ***********************************************************************
%% Added control analyses
%% ***********************************************************************




%% Control analyses for Figure 2
% Single trial traces for Ramp, PC1 and TinCon without baseline correction
load_Fig2

rows = length(dimms); colms = 2;

disp('Plotting...')
cutoff = 0.7; ls = '-'; lw = 1;
smooth_win_t = 0.025;

isNorm = 1; addN = {'', '_norm'};
sub0 = 0;

% Plot motion discrimination task
figure('Position', [200 10 600 600]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')

% Plot single trials
doMeanSub0 = 0; pltAll = 0; meanCorr = 0; doBLcorr = 0;
plot_singleTrials_highlights



%% Control analyses for Figure 6

% Mediation analysis of integral of Min_con-Min_ips
% Included in run_mediation_all2all.m
mediation_Min_spotlight


%% Control figure - average activity in ramp direction for individual sessions 

load(fullfile(saveLoc, 'sig_allSessions')) % sigSL

plotSig = {'Ramp'};
sigName = {'S^{ramp} [a.u.]'};
dimms = [63];

yLims_BLsub = {[-.6 1.6]}; % Ramp, PC1, TinC
yLims_raw = {[-.5 1.5]}; % Ramp, PC1, TinC
yTicks = {[-.5 0 .5 1]};  % Ramp, PC1, TinC

fWin = 50;
plotFiltST = gausswin(fWin, 1.5); plotFiltST = plotFiltST / sum(plotFiltST);

sessions = unique(sigSL.session)';
rows = length(sessions); colms = 2;

disp('Plotting...')
cutoff = 0.7; ls = '-'; lw = 1;
smooth_win_t = 0.025;

sub0 = 0; plotGrayLine = 1; 

% Plot motion discrimination task
% Plot trial averages
plot_resp_aligned = 1;

plot_trialAvg_decision_singleSess



