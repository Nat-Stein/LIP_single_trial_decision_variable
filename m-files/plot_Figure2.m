% plot_Figure2 of LIP single trial paper

rows = length(plotSig); colms = 4;

disp('Plotting...')
cutoff = 0.7; ls = '-'; lw = 1;
smooth_win_t = 0.025;

sub0 = 0; plotGrayLine = 1; 

% Plot motion discrimination task
figure('Position', [200 10 850 600]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')

% Plot single trials
doMeanSub0 = 0; pltAll = 0; meanCorr = 0; doBLcorr = 1;
plot_singleTrials_highlights

% Plot trial averages
plot_resp_aligned = 1;
plot_trialAvg_decision_allSess
% plot_bigSbigT_multipleS_trialAvg_decision_normS % Change away from normS


