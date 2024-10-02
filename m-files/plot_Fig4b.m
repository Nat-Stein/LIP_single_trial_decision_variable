% Plot Figure 4b

load(fullfile(saveLoc, 'decoderPerf_crossTP'))
load(fullfile(par_sess{1}.datPath, 'Plotting', 'darkRedBlue')); darkRedBlueUpper = darkRedBlue(size(darkRedBlue,1)/2+1:end, :);
% useTPTs = find(winStart >= 0.1 & winStart <= 0.5); useTPTs = useTPTs(1) : 2 : useTPTs(end);
xyLims = [0.1 0.501];  useTPTs = find(winCenter >= xyLims(1) & winCenter <= xyLims(2));
useTPTs = useTPTs(1) : 2 : useTPTs(end);
t_decoder = winCenter(useTPTs)-0.001;
figure; hold on;
plotDat = squeeze(pc_allMat(:, useTPTs, useTPTs));
imagesc(squeeze(nanmean(plotDat, 1))',[0.5,0.8]);
colormap(darkRedBlueUpper);
axis xy
yt = get(gca, 'YTick'); ytNew = yt(1:2:end-1); ytlbl = t_decoder((ytNew+1));
set(gca, 'YTick',ytNew+1, 'YTickLabel',ytlbl, 'YDir','normal')
ylabel('Test time point [s]')
xt = get(gca, 'XTick'); xtNew = xt(1:2:end-1); xtlbl = t_decoder((xtNew+1));
set(gca, 'XTick',xtNew+1, 'XTickLabel',xtlbl)
set(gca,'TickDir','out', 'FontSize', 8);
xlabel('Training time point [s]')
l1 = 16; l2 = 41; 
line([l1 l2], l1 * [1 1], 'Color', [1 1 1], 'LineStyle', '--', 'LineWidth', 1.5)
line(l1 * [1 1], [l1 l2], 'Color', [1 1 1], 'LineStyle', '--', 'LineWidth', 1.5)
colorbar
axis square
axis tight
