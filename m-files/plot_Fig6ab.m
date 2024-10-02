% plot_Fig6ab
% Plot activity of  Min neurons during motion viewing and motion
% discrimination task

if ~exist('sigSL'); load(fullfile(saveLoc, 'sig_allSessions')); end
tmp = load(fullfile(saveLoc, 'sig_view_allSessions')); sigSLv = tmp.sigSL; clear tmp
plotSig = {'MinC', 'MinI'};
sigName = {'S^{left}_{Min} [sp/s]', 'S^{right}_{Min} [sp/s]'};
rows = length(plotSig); colms = 3;
sub0 = 0; plotGrayLine = 0; plot_resp_aligned = 1;
xLimView = [0 0.4];

figure('Position', [600 400 700 600]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')
% trial averages viewing and decision - DONE
plot_trialAvg_decision_allSess

% trial averages viewing and decision
pSL = 1;
plot_trialAvg_view_allSess
