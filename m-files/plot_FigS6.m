% plot_FigS7

%% Setting up the parameters

load(fullfile(saveLoc, 'sig_allSessions'))

% For highlights plot
sigCohs = nanunique(sigSL.sig_coh); plotCohs = [6 2];
plotSig = {'WhenD', 'WhatD'};
sigName = {'S^{when} [a.u.]', 'S^{what} [a.u.]'};
yLims_BLsub = {[-.8 1.4], [-1 2]};
yTicks = {[-.5 0 .5 1], [-1 0 1 2]};

fWin = 80; fWin = 50;
plotFiltST = gausswin(fWin, 1.5); plotFiltST = plotFiltST / sum(plotFiltST);
xlims{1} = [0.2 0.5];
tBL = xlims{1}(1) + fWin/70/1000 * [-1 1]; tBLv = find(par_sess{1}.tsl >= tBL(1) & par_sess{1}.tsl < tBL(2) );

% Trials to highlight
uTRL{2} = [1 3 19 20 21 65 45];
uTRL{6} = [2 3 5 6 11 14 16 17 65];
% Trials not to display
exclTRL{2} = [17 79];
exclTRL{6} = [4 24 52];
pltDash = {'-', [], ':'};

numTr = 80;
choiceCol = 1; choiceCols = {[0 0 0], [201 178 4]/255};
cohText{2} = {'25.6% contra'};
cohText{6} = {'0%'};


%% Plotting

rows = length(plotSig); colms = 4;

disp('Plotting...')
cutoff = 0.7; ls = '-'; lw = 1;
smooth_win_t = 0.025;

sub0 = 0; plotGrayLine = 1; 

% Plot motion discrimination task
figure('Position', [200 10 850 450]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')

% Plot single trials
doMeanSub0 = 0; pltAll = 0; meanCorr = 0; doBLcorr = 1;
plot_singleTrials_highlights

% Plot trial averages
plot_resp_aligned = 1;
plot_trialAvg_decision_allSess
% plot_bigSbigT_multipleS_trialAvg_decision_normS % Change away from normS


