% Plot stimulus-locked mean traces during viewing task for all coherences 
% - considering nans in each trial
% meanTracesSLview_bigSbigT

filter = ones(size(sigSLv.sig_coh))==1;
conditions = sigSLv.sig_coh;
allCond = nanunique(sigSL.sig_coh);
colores = NS_colorsBY(length(allCond));

clear Son
dt = sigSL.t(2)-sigSL.t(1);
min_num_bins = 50; % read from elsewhere...
[Son,ti_on] = bigSbigT(sl,dt,dt,smooth_win_t,conditions,filter,min_num_bins);
ton = sigSL.t(ti_on);

sCoh = nanunique(conditions)';

if ~exist('sub0'); sub0 = 0; end

coh0 = find(sCoh == 0);
maxY = 0; minY = 0;
for c = 1 : length(sCoh)
    
    if sub0 == 0
        son = Son(c, :);
    else
        son = Son(c, :) - Son(coh0, :);
    end
    coh = sCoh(c);
    trls = find(conditions == coh &...
        filter == 1);
    
    nans_on  = nanmean(isnan(sl(trls,ti_on)));
    tind_on  = nans_on<cutoff;
    
    I1 = tind_on & ~isnan(son);
    t_plt = ton(I1);
    s_plt = son(I1);
    
    col = find(allCond == coh);
    plot(t_plt, s_plt,...
        'Color', colores(col, :), 'LineStyle', ls, 'LineWidth', lw)
    
    maxY = max(maxY, max(s_plt));
    minY = min(minY, min(s_plt));
    
end




