% plot_FigS3

load(fullfile(saveLoc, 'sig_allSessions'))

dimms = [63 20 1]; 
plotSig = {'Ramp', 'PC1', 'TinC'};
sigName = {'S^{ramp} [a.u.]', 'S^{PC1} [a.u.]', 'S^{con}_{Tin} [sp/s]'};


yTicks = {[0 .5 1], [0 .5 1], [0 10 20 30]}; 
% cutoff = 0.7; ls = '-'; lw = 1;
% smooth_win_t = 0.025;
rows = length(dimms); colms = 2; 

saveFig = 0;

sub0 = 0; plotGrayLine = 0; 
plot_resp_aligned = 1;

nQuant = 4; 
qRT = quantile(sigSL.rt, nQuant); qRT = [0 qRT max(sigSL.rt)];
rt_quant = nan(size(sigSL.rt)); 
for qu = 1 : length(qRT)-1
    qq = find((sigSL.rt >= qRT(qu) & sigSL.rt <= qRT(qu+1)) & sigSL.choice == 1);
    rt_quant(qq) = qu; 
    qq = find((sigSL.rt >= qRT(qu) & sigSL.rt <= qRT(qu+1)) & sigSL.choice == 0);
    rt_quant(qq) = qu+nQuant+1; 
end

% Plot motion discrimination task split by RT percentile
figure('Position', [200 10 300 600]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')
i=0;
plot_bigSbigT_multipleSignals_decisionRT


















