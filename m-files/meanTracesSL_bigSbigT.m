% Plot stimulus-locked mean traces for all coherences - considering nans in
% each trial
% meanTracesSL_bigSbigT

if exist('sigSL') == 0 && exist('normS') >= 1
    sigSL.sig_coh = normS.sigSL.sig_coh;
    sigSL.t = normS.sigSL.t;
    
elseif exist('sigSL') == 0 && exist('normS') == 0
    error('No signal structure loaded')
end

filter = ones(size(sigSL.sig_coh))==1;
if exist('singleSess')
    if singleSess == 1
        filter = sigSL.session == sess;
    end
end
conditions = sigSL.sig_coh; % conditions = adummyvar(conditions);
colores = NS_colorsBY(length(nanunique(conditions)));

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
    
    plot(t_plt, s_plt,...
        'Color', colores(c, :), 'LineStyle', ls, 'LineWidth', lw)
    
    maxY = max(maxY, max(s_plt));
    minY = min(minY, min(s_plt));
    
end




