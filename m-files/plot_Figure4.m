% plot_Figure4: The population signal predictive of choice and RT is approximately one-dimensional.


% a) Decoder performance over time
% variables computed in: calc_whatDecoder_fixed_time
plot_Fig4a
% Is missing simulation!

% b) Decoder percormance across time points
plot_Fig4b

% c) Trial averages for when direction
if ~exist('sigSL'); load(fullfile(saveLoc, 'sig_allSessions')); end
plotSig = {'WhenD'};
sigName = {'S^{when} [a.u.]'};

rows = length(plotSig); colms = 2;
sub0 = 0; plotGrayLine = 0; 
plot_resp_aligned = 1;

figure('Position', [600 400 700 300]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')
% Plot trial averages
plot_trialAvg_decision_allSess

% d) Cosine similarity & Within-trial correlation
useD = [3 4 1 5 6]; 
dimLabelsCS = {'ramp', 'PC1', 'Tin^{con}_{in}', 'What-decoder', 'When-decoder'};
% Compute cosine similarities between weight vectors
calc_weightCS
% Compute within-trial correlation between signals
calc_wiTrialCorr
% Plot mean CS and within-trial correlation
%plot_Fig4de 
plot_Fig4de_az