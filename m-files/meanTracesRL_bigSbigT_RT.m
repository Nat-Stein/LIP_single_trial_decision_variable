% Plot stimulus-locked mean traces for all coherences - considering nans in
% each trial
% meanTracesRL_bigSbigT


conditions = rt_quant;
colores = NS_colorsBY(length(nanunique(conditions)));

dt = sigRL.t(2)-sigRL.t(1);
min_num_bins = 50; % read from elsewhere...
[Son,ti_on] = bigSbigT(rl,dt,dt,smooth_win_t,conditions,filter,min_num_bins);
ton = sigRL.t(ti_on);

quant = nanunique(conditions)';


for c = 1 : length(quant)
    
    son = Son(c, :);
    qnt = quant(c);
    trls = find(conditions == qnt &...
        filter == 1);
    
    nans_on  = nanmean(isnan(rl(trls,ti_on)));
    tind_on  = nans_on<cutoff;
    
    I1 = tind_on & ~isnan(son);
    t_plt = ton(I1);
    s_plt = son(I1);
    
    plot(t_plt, s_plt,...
        'Color', colores(c, :), 'LineStyle', ls, 'LineWidth', lw)
    
    maxY = max(maxY, max(s_plt));
    minY = min(minY, min(s_plt));
    
end




