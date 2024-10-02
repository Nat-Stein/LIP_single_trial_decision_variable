% plot_FigS4 of LIP single trial paper

%% load signals for plotting single trials and trial averages

load(fullfile(saveLoc, 'sig_allSessions'))

plotSig = {'Ramp', 'PC1', 'TinC'};
sigName = {'S^{ramp} [a.u.]', 'S^{PC1} [a.u.]', 'S^{con}_{Tin} [sp/s]'};


sigCohs = nanunique(sigSL.sig_coh); 

fWin = 80; fWin = 50;
plotFiltST = gausswin(fWin, 1.5); plotFiltST = plotFiltST / sum(plotFiltST);
xlims{1} = [0.2 0.5];
tBL = xlims{1}(1) + fWin/70/1000 * [-1 1]; tBLv = find(par_sess{1}.tsl >= tBL(1) & par_sess{1}.tsl < tBL(2) );

%% Plot 

rows = length(sigName); colms = 1;

disp('Plotting...')
cutoff = 0.7; ls = '-'; lw = 1;
smooth_win_t = 0.025;

sub0 = 1; plotGrayLine = 0; 

% Plot motion discrimination task
figure('Position', [200 10 300 600]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')

doMeanSub0 = 0; pltAll = 0; meanCorr = 0; doBLcorr = 1;
plot_resp_aligned = 0; 

% Plot trial averages
plot_trialAvg_decision_allSess


