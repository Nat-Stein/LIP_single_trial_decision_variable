% Plot stimulus-locked mean traces for all coherences - considering nans in
% each trial
% meanTracesRL_bigSbigT

if exist('sigRL') == 0 && exist('normS') >= 1
    sigRL.sig_coh = normS.sigRL.sig_coh;
    sigRL.t = normS.sigRL.t;
    sigRL.correct = normS.sigRL.correct;
    sigRL.session = normS.sigRL.session;
    sigRL.choice = normS.sigRL.choice;
elseif exist('sigRL') == 0 && exist('normS') == 0
    error('No signal structure loaded')
end


filter = sigRL.correct==1;
if exist('singleSess')
    if singleSess == 1
        filter = sigSL.session == sess & sigRL.correct==1;
    end
end
conditions = sigRL.sig_coh;
I = sigRL.sig_coh==0 & ((sigRL.choice==1 & sigRL.correct==1) | (sigRL.choice==0 & sigRL.correct==0));
conditions(I) = 0.0001;
I = sigRL.sig_coh==0 & ((sigRL.choice==0 & sigRL.correct==1) | (sigRL.choice==1 & sigRL.correct==0));
conditions(I) = -0.0001;
 % conditions = adummyvar(conditions);

colores = NS_colorsBY(length(nanunique(conditions)));

dt = sigRL.t(2)-sigRL.t(1);
min_num_bins = 50; % read from elsewhere...
[Son,ti_on] = bigSbigT(rl,dt,dt,smooth_win_t,conditions,filter,min_num_bins);
ton = sigRL.t(ti_on);

sCoh = nanunique(conditions)';


for c = 1 : length(sCoh)
    
    son = Son(c, :);
    coh = sCoh(c);
    trls = find(conditions == coh &...
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




