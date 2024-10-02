% load_Fig2
% load signals for plotting single trials and trial averages

load(fullfile(saveLoc, 'sig_allSessions'))

plotSig = {'Ramp', 'PC1', 'TinC'};
sigName = {'S^{ramp} [a.u.]', 'S^{PC1} [a.u.]', 'S^{con}_{Tin} [sp/s]'};
dimms = [63 20 61];

% From plot_singleTrials_highlights 


% For highlights plot
sigCohs = nanunique(sigSL.sig_coh); plotCohs = [6 2];
sigName = {'S^{ramp} [a.u.]', 'S^{PC1} [a.u.]', 'S^{con}_{Tin} [sp/s]'};
% yLims = {[0.3 2], [.3 2.2], [5 75]}; % Ramp, PC1, TinC
yLims_BLsub = {[-.6 1.4], [-.6 1.4], [-25 40]}; % Ramp, PC1, TinC
yLims_raw = {[-.5 1.3], [-.5 1.3], [0 50]}; % Ramp, PC1, TinC
yTicks = {[-.5 0 .5 1], [-.5 0 .5 1], [-20 0 20 40]};  % Ramp, PC1, TinC

fWin = 50;
plotFiltST = gausswin(fWin, 1.5); plotFiltST = plotFiltST / sum(plotFiltST);
xlims{1} = [0.2 0.5];
% tBL = xlims{1}(1) + fWin/70/1000 * [-1 1]; tBLv = find(par_sess{1}.tsl >= tBL(1) & par_sess{1}.tsl < tBL(2) );
tBLv = find(round(par_sess{1}.tsl, 3) == round(xlims{1}(1), 3)); 

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