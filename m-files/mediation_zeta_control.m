% mediation_zeta_control

load(fullfile(saveLoc, 'sig_allSessions'),...
    'sigSL')
S = sigSL;
% Analysis parameters
max_coh = 0.1;
minRT = 0.67;
maxRT = 2;
start_t = 0.39;
end_t = 0.5;

session = S.session;
times = S.t;
dt = times(2) - times(1);
choice = S.choice;
RT = S.rt;
coh = S.sig_coh;

str =  {'PC1',...
    'TinC',...
    'ramp'};

% Run mediation only on one time point!


plotFigs = 0;

% smooth in window
for i = 1 : length(str)
    sm = round(0.05/dt);
    h = ones(sm,1)/sm;
    S.(str{i}) = conv2(1, h, S.(str{i}), 'same');
end

%% time subset
tind = findclose(times, -0.1:0.01:0.75);
for i=1:length(str)
    aux = S.(str{i});
    S.(str{i}) = aux(:,tind);
end
times = times(tind);


%% Run mediation itself

warning('off','all')

flags.norm_to_se = 1;
nsessions =  8;
clear lev_ctrl
count = 0; 
for iMethod=1:length(str)
    sProj = S.(str{iMethod});
    disp(['Method ' num2str(iMethod) '/' num2str(length(str))]);
    for iMediator=1:length(str)
        
        if iMediator~=iMethod
            
            tind = findclose(times, 0.55);
            sProj_mediator = sProj(:,tind);
            sOther_mediator = S.(str{iMediator})(:,tind);
            
            count = count + 1;
            
            %% per session and average
            do_plot = 0;
            clear out out_other
            clear thisRT thisRT_part thisRT_part_other
            clear thisCho thisCho_part thisCho_part_other
            clear thisRT_perm thisRT_part_perm thisRT_part_other_perm
            clear thisCho_perm thisCho_part_perm thisCho_part_other_perm
            i = 0;
            for j=1:nsessions
                I = session==j & abs(coh)<max_coh; % only low coh
                if ~[all(isnan(to_vec(sProj(I,:)))) || all(isnan(sProj_mediator(I))) || all(isnan(sOther_mediator(I)))]
                    i = i+1;
                    disp(['Session ' num2str(i)])
                    
                    sProj_orig = sProj(I,:);
                    sProj_mediator_orig = sProj_mediator(I);
                    sOther_mediator_orig = sOther_mediator(I);
                    choice_orig = choice(I);
                    RT_orig = RT(I);
                    coh_orig = coh(I);
                    useCoh = unique(coh_orig)';
                    numTR = length(coh_orig);
                    orig_trl = 1 : numTR;
                    
                    
                    for perm = 1 : 1000
                        
                        k = datasample(orig_trl, numTR);
                        [~,out_samp(perm)] = corr_with_RT_choice(choice_orig(k)==0, RT_orig(k), coh_orig(k), times, ...
                            sProj_orig(k,:), sProj_mediator_orig(k), minRT,maxRT, start_t, end_t, do_plot,flags);
                        [~,out_other_samp(perm)] = corr_with_RT_choice(choice_orig(k)==0, RT_orig(k), coh_orig(k), times, ...
                            sProj_orig(k,:), sOther_mediator_orig(k), minRT,maxRT, start_t, end_t,do_plot,flags);
                        
                        coh_k = coh_orig(k)';
                        m = nan(size(k));
                        for cc = useCoh
                            cohtrl = find(coh_k == cc);
                            m(cohtrl) =  k(cohtrl(randperm(length(cohtrl))));
                        end
                        m = k(randperm(numTR));
                        [~,out_perm(perm)] = corr_with_RT_choice(choice_orig(k)==0, RT_orig(k), coh_orig(k), times, ...
                            sProj_orig(k,:), sProj_mediator_orig(m), minRT,maxRT, start_t, end_t, do_plot,flags);
                        [~,out_other_perm(perm)] = corr_with_RT_choice(choice_orig(k)==0, RT_orig(k), coh_orig(k), times, ...
                            sProj_orig(k,:), sOther_mediator_orig(m), minRT,maxRT, start_t, end_t,do_plot,flags);
                        
                        thisRT(:, perm, i) = out_samp(perm).rho_RT;
                        thisRT_part(:, perm, i) = out_samp(perm).rho_RT_partial;
                        thisRT_part_other(:, perm, i) = out_other_samp(perm).rho_RT_partial;
                        
                        thisCho(:, perm, i) = out_samp(perm).rho_choice;
                        thisCho_part(:, perm, i) = out_samp(perm).rho_choice_partial;
                        thisCho_part_other(:, perm, i) = out_other_samp(perm).rho_choice_partial;
                        
                        thisRT_perm(:, perm, i) = out_perm(perm).rho_RT;
                        thisRT_part_perm(:, perm, i) = out_perm(perm).rho_RT_partial;
                        thisRT_part_other_perm(:, perm, i) = out_other_perm(perm).rho_RT_partial;
                        
                        thisCho_perm(:, perm, i) = out_perm(perm).rho_choice;
                        thisCho_part_perm(:, perm, i) = out_perm(perm).rho_choice_partial;
                        thisCho_part_other_perm(:, perm, i) = out_other_perm(perm).rho_choice_partial;
                        
                    end
                end
            end
            
            %% save
            
            lev_ctrl(count).signal = str{iMethod};
            lev_ctrl(count).mediator = str{iMediator};
            
            lev_ctrl(count).t = out_samp(1).tt;
            
            lev_ctrl(count).RT.sig.rho = thisRT;
            lev_ctrl(count).RT.sig.rho_part = thisRT_part;
            lev_ctrl(count).RT.sig.rho_part_other = thisRT_part_other;
            
            lev_ctrl(count).choice.sig.rho = thisCho;
            lev_ctrl(count).choice.sig.rho_part = thisCho_part;
            lev_ctrl(count).choice.sig.rho_part_other = thisCho_part_other;
            
            lev_ctrl(count).RT.perm.rho = thisRT_perm;
            lev_ctrl(count).RT.perm.rho_part = thisRT_part_perm;
            lev_ctrl(count).RT.perm.rho_part_other = thisRT_part_other_perm;
            
            lev_ctrl(count).choice.perm.rho = thisCho_perm;
            lev_ctrl(count).choice.perm.rho_part = thisCho_part_perm;
            lev_ctrl(count).choice.perm.rho_part_other = thisCho_part_other_perm;
            
            
            save(fullfile(saveLoc, 'lev_ctrl_zeta2'), 'lev_ctrl', 'start_t', 'end_t')

        end
    end
    
end

save(fullfile(saveLoc, 'lev_ctrl_zeta2'), 'lev_ctrl', 'start_t', 'end_t')


% Correct to shuffle within coherence!!!

%% Compute zetas
% cts = [1 5 8]; ct_tin = 5;

cts = [1 3 6]; ct_tin = 3;

tpt = find(round(lev_ctrl(1).t, 3) == 0.4);
med = [];

for count = cts
    disp([num2str(count) ' ' lev_ctrl(count).signal ' mediated by ' lev_ctrl(count).mediator])
    
    % ---------------------------------------
    % RT - data
    raw = squeeze(lev_ctrl(count).RT.sig.rho(tpt, :, :));
    
    % by self
    partial = squeeze(lev_ctrl(count).RT.sig.rho_part(tpt, :, :));
    partial(partial>0) = 0;
    medRT = 1 - partial.^2./raw.^2;
    medRT(medRT<0) = 0;
    med(count).RT.sig.self = medRT;
    
    % by Tin
    partial_tin = squeeze(lev_ctrl(count).RT.sig.rho_part_other(tpt, :, :));
    partial_tin(partial_tin>0) = 0;
    medRT = 1 - partial_tin.^2./raw.^2;
    medRT(medRT<0) = 0;
    med(count).RT.sig.tin = medRT;
    
    % ---------------------------------------
    % RT - perm
    raw = squeeze(lev_ctrl(count).RT.perm.rho(tpt, :, :));
    
    % by self
    partial = squeeze(lev_ctrl(count).RT.perm.rho_part(tpt, :, :));
    partial(partial>0) = 0;
    medRT = 1 - partial.^2./raw.^2;
    medRT(medRT<0) = 0;
    med(count).RT.perm.self = medRT;
    
    % by Tin
    partial_tin = squeeze(lev_ctrl(count).RT.perm.rho_part_other(tpt, :, :));
    partial_tin(partial_tin>0) = 0;
    medRT = 1 - partial_tin.^2./raw.^2;
    medRT(medRT<0) = 0;
    med(count).RT.perm.tin = medRT;
    
    % ---------------------------------------
    
     % ---------------------------------------
    % choice - data
    raw = squeeze(lev_ctrl(count).choice.sig.rho(tpt, :, :));
    
    % by self
    partial = squeeze(lev_ctrl(count).choice.sig.rho_part(tpt, :, :));
    medCho = nan(size(partial));
    medCho = 1 - partial ./ raw;
    medCho(partial<=0 & raw > 0) = 1;
    medCho(raw<=0) = nan;
    medCho(medCho<0) = 0;
    med(count).choice.sig.self = medCho;
    
    % by Tin
    partial_tin = squeeze(lev_ctrl(count).choice.sig.rho_part_other(tpt, :, :));
    medCho = nan(size(partial_tin));
    medCho = 1 - partial_tin ./ raw;
    medCho(partial_tin<=0 & raw > 0) = 1;
    medCho(raw<=0) = nan;
    medCho(medCho<0) = 0;
    med(count).choice.sig.tin = medCho;
    
    % ---------------------------------------
    % choice - perm
    raw = squeeze(lev_ctrl(count).choice.perm.rho(tpt, :, :));
    
    % by self
    partial = squeeze(lev_ctrl(count).choice.perm.rho_part(tpt, :, :));
    medCho = nan(size(partial));
    medCho = 1 - partial ./ raw;
    medCho(partial<=0 & raw > 0) = 1;
    medCho(raw<=0) = nan;
    medCho(medCho<0) = 0;
    med(count).choice.perm.self = medCho;
    
    % by Tin
    partial_tin = squeeze(lev_ctrl(count).choice.perm.rho_part_other(tpt, :, :));
    medCho = nan(size(partial_tin));
    medCho = 1 - partial_tin ./ raw;
    medCho(partial_tin<=0 & raw > 0) = 1;
    medCho(raw<=0) = nan;
    medCho(medCho<0) = 0;
    med(count).choice.perm.tin = medCho;
    
    % ---------------------------------------
    % Stats within session
    for sess = 1 : 8
        
        % Wilcoxon rank sum test
        
        % RT - self
        x = med(count).RT.sig.self(:, sess);
        y = med(count).RT.perm.self(:, sess);
        [pval,h,stats] = ranksum(x,y);
        pvals.RT(count, 1, sess) = pval;
        zvals.RT(count, 1, sess) = stats.zval;
        
        % RT - tin
        x = med(count).RT.sig.tin(:, sess);
        y = med(count).RT.perm.tin(:, sess);
        [pval,h,stats] = ranksum(x,y);
        pvals.RT(count, 2, sess) = pval;
        zvals.RT(count, 2, sess) = stats.zval;
        
        % Choice - self
        x = med(count).choice.sig.self(:, sess);
        y = med(count).choice.perm.self(:, sess);
        [pval,h,stats] = ranksum(x,y);
        pvals.choice(count, 1, sess) = pval;
        zvals.choice(count, 1, sess) = stats.zval;
        % Choice - tin
        x = med(count).choice.sig.tin(:, sess);
        y = med(count).choice.perm.tin(:, sess);
        [pval,h,stats] = ranksum(x,y);
        pvals.choice(count, 2, sess) = pval;
        zvals.choice(count, 2, sess) = stats.zval;
    end
    disp(['RT by self, max p = ' num2str(max(max(squeeze(pvals.RT(count, 1, :)))))])
    disp(['Choice by self, max p = ' num2str(max(max(squeeze(pvals.choice(count, 1, :)))))])
    if count ~= ct_tin
        disp(['RT by Tin, max p = ' num2str(max(max(squeeze(pvals.RT(count, 2, :)))))])
        disp(['Choice by Tin, max p = ' num2str(max(max(squeeze(pvals.choice(count, 2, :)))))])
    end
    disp('---------------------')
    % ---------------------------------------
    % Stats across sessions
    x = mean(med(count).RT.sig.self, 2);
    y = mean(med(count).RT.perm.self, 2);
    [pval,h,stats] = ranksum(x,y);
    pvals.RTm(count, 1) = pval;
    zvals.RTm(count, 1) = stats.zval;
    
    % RT - tin
    x = mean(med(count).RT.sig.tin, 2);
    y = mean(med(count).RT.perm.tin, 2);
    [pval,h,stats] = ranksum(x,y);
    pvals.RTm(count, 2) = pval;
    zvals.RTm(count, 2) = stats.zval;
    
    % Choice - self
    x = mean(med(count).choice.sig.self, 2);
    y = mean(med(count).choice.perm.self, 2);
    [pval,h,stats] = ranksum(x,y);
    pvals.choicem(count, 1) = pval;
    zvals.choicem(count, 1) = stats.zval;
    % Choice - tin
    x = mean(med(count).choice.sig.tin, 2);
    y = mean(med(count).choice.perm.tin, 2);
    [pval,h,stats] = ranksum(x,y);
    pvals.choicem(count, 2) = pval;
    zvals.choicem(count, 2) = stats.zval;
    
    disp(['RT mean by self, max p = ' num2str(max(squeeze(pvals.RTm(count, 1))))])
    disp(['Choice mean by self, max p = ' num2str(max(squeeze(pvals.choicem(count, 1))))])
    if count ~= ct_tin
        disp(['RT mean by Tin, max p = ' num2str(max(squeeze(pvals.RTm(count, 2))))])
        disp(['Choice mean by Tin, max p = ' num2str(max(squeeze(pvals.choicem(count, 2))))])
    end
    
    disp('---------------------')
    disp('---------------------')
end

all_self_RT = squeeze(pvals.RT(cts, 1, :)); median(all_self_RT(:))
all_self_choice = squeeze(pvals.choice(cts, 1, :)); median(all_self_choice(:))

allRT = all_self_RT(:)';
allChoice = all_self_choice(:)';

median([all_self_choice(:)' all_self_RT(:)'])


cts_tin = [1 8];
all_tin_RT = squeeze(pvals.RT(cts_tin, 2, :)); median(all_tin_RT(:))
all_tin_choice = squeeze(pvals.choice(cts_tin, 2, :)); median(all_tin_choice(:))

allRT = all_tin_RT(:)';
allChoice = all_tin_choice(:)';

median([all_tin_choice(:)' all_tin_RT(:)'])

% % 1 PC1 mediated by TinC
% RT by self, max p = 1.3414e-62
% Choice by self, max p = 3.1416e-106
% RT by Tin, max p = 1.5059e-12
% Choice by Tin, max p = 1.1325e-28
% ---------------------
% RT mean by self, max p = 0
% Choice mean by self, max p = 0
% RT mean by Tin, max p = 0
% Choice mean by Tin, max p = 0
% ---------------------
% ---------------------
% % 5 TinC mediated by PC1
% RT by self, max p = 0.022617
% Choice by self, max p = 0.7343
% ---------------------
% RT mean by self, max p = 0
% Choice mean by self, max p = 6.7921e-307
% ---------------------
% ---------------------
% % 8 ramp mediated by TinC
% RT by self, max p = 3.2157e-182
% Choice by self, max p = 1.636e-223
% RT by Tin, max p = 0.036706
% Choice by Tin, max p = 7.6914e-34
% ---------------------
% RT mean by self, max p = 0
% Choice mean by self, max p = 5.5169e-258
% RT mean by Tin, max p = 0
% Choice mean by Tin, max p = 5.0105e-258
% ---------------------
% ---------------------

% Rerun with permutation within coherence
% 1 PC1 mediated by TinC
% RT by self, max p = 9.4504e-68
% Choice by self, max p = 6.4884e-94
% RT by Tin, max p = 4.2162e-17
% Choice by Tin, max p = 8.6409e-36
% ---------------------
% RT mean by self, max p = 0
% Choice mean by self, max p = 0
% RT mean by Tin, max p = 0
% Choice mean by Tin, max p = 0
% ---------------------
% ---------------------
% 3 TinC mediated by PC1
% RT by self, max p = 0.00021795
% Choice by self, max p = 0.90508
% ---------------------
% RT mean by self, max p = 0
% Choice mean by self, max p = 1.2628e-306
% ---------------------
% ---------------------
% 6 ramp mediated by TinC
% RT by self, max p = 2.9401e-169
% Choice by self, max p = 6.1959e-244
% RT by Tin, max p = 0.047021
% Choice by Tin, max p = 1.0297e-42
% ---------------------
% RT mean by self, max p = 0
% Choice mean by self, max p = 4.1308e-261
% RT mean by Tin, max p = 0
% Choice mean by Tin, max p = 4.0378e-261
% ---------------------
% ---------------------

% Within individual sessions, choice
count = 3;
hb = [0 : 0.01 : 1];
figure; hold on
for sess = 1 : 8
    subplot(2, 4, sess); hold on
    a = hist(med(count).choice.sig.self(:, sess), hb);
    b = hist(med(count).choice.perm.self(:, sess), hb);
    plot(hb, a)
    plot(hb, b)
    ylim([0 50])
    title(['S ' num2str(sess) ', p=' num2str(pvals.choice(count, 1, sess))])
    ylabel('# permutations')
    xlabel('zeta')
end
legend('signal', 'control')


%% Compare zeta real to zeta control distrbutions

nPerm = 1000;

for count = cts
    disp([num2str(count) ' ' lev_ctrl(count).signal ' mediated by ' lev_ctrl(count).mediator])
    
    % ---------------------------------------
    % RT - self
    data1 = nanmean(med(count).RT.sig.self);
    data2 = nanmean(med(count).RT.perm.self);
    SE1 = nanstd(med(count).RT.sig.self)/sqrt(nPerm);
    SE2 = nanstd(med(count).RT.perm.self/sqrt(nPerm));
    
    [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
    
    pairedT_zeta(count).RT.tStat = tStat;
    pairedT_zeta(count).RT.pValue = pValue;
    
    disp(['RT_self: t=' num2str(tStat) ', p=' num2str(pValue)])
    
    % ---------------------------------------
    % Choice - self
    data1 = nanmean(med(count).choice.sig.self);
    data2 = nanmean(med(count).choice.perm.self);
    SE1 = nanstd(med(count).choice.sig.self)/sqrt(nPerm);
    SE2 = nanstd(med(count).choice.perm.self)/sqrt(nPerm);
    
    [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
    
    pairedT_zeta(count).choice.tStat = tStat;
    pairedT_zeta(count).choice.pValue = pValue;
    
    disp(['Choice_self: t=' num2str(tStat) ', p=' num2str(pValue)])
    
    if count ~= ct_tin
        % ---------------------------------------
        % RT - Tin
        data1 = nanmean(med(count).RT.sig.tin);
        data2 = nanmean(med(count).RT.perm.tin);
        SE1 = nanstd(med(count).RT.sig.tin)/sqrt(nPerm);
        SE2 = nanstd(med(count).RT.perm.tin)/sqrt(nPerm);
        
        [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
        
        pairedT_zeta(count).RT_tin.tStat = tStat;
        pairedT_zeta(count).RT_tin.pValue = pValue;
        
        disp(['RT_Tin: t=' num2str(tStat) ', p=' num2str(pValue)])
        
        % ---------------------------------------
        % Choice - Tin
        data1 = nanmean(med(count).choice.sig.tin);
        data2 = nanmean(med(count).choice.perm.tin);
        SE1 = nanstd(med(count).choice.sig.tin)/sqrt(nPerm);
        SE2 = nanstd(med(count).choice.perm.tin)/sqrt(nPerm);
        
        [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
        
        pairedT_zeta(count).choice_tin.tStat = tStat;
        pairedT_zeta(count).choice_tin.pValue = pValue;
        
        disp(['Choice_Tin: t=' num2str(tStat) ', p=' num2str(pValue)])
        
    end
    
    disp('-----')
end

% Using SE
% 1 PC1 mediated by TinC
% RT_self: t=-4317.3372, p=9.4463e-24
% Choice_self: t=-840.4379, p=8.9171e-19
% RT_Tin: t=-2116.3909, p=1.3887e-21
% Choice_Tin: t=-536.4806, p=2.0647e-17
% -----
% 3 TinC mediated by PC1
% RT_self: t=-2378.0673, p=6.1405e-22
% Choice_self: t=-1619.4182, p=9.0419e-21
% -----
% 6 ramp mediated by TinC
% RT_self: t=-5720.6593, p=1.3172e-24
% Choice_self: t=-2050.5185, p=1.7327e-21
% RT_Tin: t=-1251.1454, p=5.5032e-20
% Choice_Tin: t=-433.3177, p=9.2062e-17
% 
% % Using std
% 1 PC1 mediated by TinC
% RT_self: t=-136.5262, p=2.9837e-13
% Choice_self: t=-26.577, p=2.7347e-08
% RT_Tin: t=-66.9262, p=4.3701e-11
% Choice_Tin: t=-16.965, p=6.0593e-07
% -----
% 3 TinC mediated by PC1
% RT_self: t=-75.2011, p=1.9343e-11
% Choice_self: t=-51.2105, p=2.8357e-10
% -----
% 6 ramp mediated by TinC
% RT_self: t=-180.9031, p=4.1626e-14
% Choice_self: t=-64.8431, p=5.4511e-11
% RT_Tin: t=-39.5647, p=1.7163e-09
% Choice_Tin: t=-13.7027, p=2.5981e-06



%% Compare raw vs partial distributions

for count = cts
    disp([num2str(count) ' ' lev_ctrl(count).signal ' mediated by ' lev_ctrl(count).mediator])
    
    % ---------------------------------------
    % RT - self
    raw = squeeze(lev_ctrl(count).RT.sig.rho(tpt, :, :));
    part = squeeze(lev_ctrl(count).RT.sig.rho_part(tpt, :, :));
    partT = squeeze(lev_ctrl(count).RT.sig.rho_part_other(tpt, :, :));
    
    data1 = nanmean(raw);
    data2 = nanmean(part);
    SE1 = nanstd(raw)/sqrt(nPerm);
    SE2 = nanstd(part)/sqrt(nPerm);
    
    [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
    
    pairedT_part(count).RT.tStat = tStat;
    pairedT_part(count).RT.pValue = pValue;
    
    disp(['RT_self: t=' num2str(tStat) ', p=' num2str(pValue)])
    
    if count ~= ct_tin
        % ---------------------------------------
        % RT - Tin
        
        data1 = nanmean(raw);
        data2 = nanmean(partT);
        SE1 = nanstd(raw)/sqrt(nPerm);
        SE2 = nanstd(part)/sqrt(nPerm);
        
        [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
        
        pairedT_partTin(count).RT.tStat = tStat;
        pairedT_partTin(count).RT.pValue = pValue;
        
        disp(['RT_Tin: t=' num2str(tStat) ', p=' num2str(pValue)])
    end
    
    % ---------------------------------------
    % RT - self
    raw = squeeze(lev_ctrl(count).choice.sig.rho(tpt, :, :));
    part = squeeze(lev_ctrl(count).choice.sig.rho_part(tpt, :, :));
    partT = squeeze(lev_ctrl(count).choice.sig.rho_part_other(tpt, :, :));
    
    data1 = nanmean(raw);
    data2 = nanmean(part);
    SE1 = nanstd(raw)/sqrt(nPerm);
    SE2 = nanstd(part)/sqrt(nPerm);
    
    [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
    
    pairedT_part(count).choice.tStat = tStat;
    pairedT_part(count).choice.pValue = pValue;
    
    disp(['choice_self: t=' num2str(tStat) ', p=' num2str(pValue)])
    
    if count ~= ct_tin
        % ---------------------------------------
        % RT - Tin
        
        data1 = nanmean(raw);
        data2 = nanmean(partT);
        SE1 = nanstd(raw)/sqrt(nPerm);
        SE2 = nanstd(part)/sqrt(nPerm);
        
        [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2);
        
        pairedT_partTin(count).choice.tStat = tStat;
        pairedT_partTin(count).choice.pValue = pValue;
        
        disp(['choice_Tin: t=' num2str(tStat) ', p=' num2str(pValue)])
    end
    
    
    disp('-----')
end 

% Use SE
% 1 PC1 mediated by TinC
% RT_self: t=2836.2439, p=1.7888e-22
% RT_Tin: t=1919.8291, p=2.7475e-21
% choice_self: t=-72.1856, p=2.5752e-11
% choice_Tin: t=-53.8508, p=1.9961e-10
% -----
% 3 TinC mediated by PC1
% RT_self: t=2278.0386, p=8.2956e-22
% choice_self: t=-86.3145, p=7.3777e-12
% -----
% 6 ramp mediated by TinC
% RT_self: t=4447.3044, p=7.6753e-24
% RT_Tin: t=1971.1884, p=2.2839e-21
% choice_self: t=-129.7245, p=4.2664e-13
% choice_Tin: t=-54.2105, p=1.9054e-10
% 
% 
% 
% 1 PC1 mediated by TinC
% RT_self: t=89.6899, p=5.6415e-12
% RT_Tin: t=60.7103, p=8.6372e-11
% choice_self: t=-2.2827, p=0.056411
% choice_Tin: t=-1.7029, p=0.13237
% -----
% 3 TinC mediated by PC1
% RT_self: t=72.0379, p=2.6123e-11
% choice_self: t=-2.7295, p=0.02936
% -----
% 6 ramp mediated by TinC
% RT_self: t=140.6361, p=2.4245e-13
% RT_Tin: t=62.3345, p=7.182e-11
% choice_self: t=-4.1022, p=0.0045601
% choice_Tin: t=-1.7143, p=0.1302


%% Raw measure and partial with SE measured directly

medN = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se_MinCI_BL'));

cnts = [25];
tpt_N = find(medN.leverage(count).t == 0.4);

for count = cnts
    disp([num2str(count) ' ' medN.leverage(count).signal ' mediated by ' medN.leverage(count).mediator])
    for sess = 1 : 8
        
        medN.leverage(count).per_session.self_mediated(sess).rho_choice(tpt_N)
        medN.leverage(count).per_session.self_mediated(sess).norm_choice(tpt_N)
        
    end
end









%%




figure; hold on
for sess = 1 : 8
    subplot(2, 4, sess); hold on
    for perm = 1 : 3
        plot(out_samp(perm).tt, thisRT(:, perm, sess), 'k')
        plot(out_samp(perm).tt, thisRT_part(:, perm, sess), 'r')
        plot(out_samp(perm).tt, thisRT_part_other(:, perm, sess), 'g')
    end
end


figure; hold on
for sess = 1 : 8
    subplot(2, 4, sess); hold on
    for perm = 1 : 3
        plot(out_samp(perm).tt, thisRT_perm(:, perm, sess), 'k')
        plot(out_samp(perm).tt, thisRT_part_perm(:, perm, sess), 'r')
        plot(out_samp(perm).tt, thisRT_part_other_perm(:, perm, sess), 'g')
    end
end

figure; hold on
for sess = 1 : 8
    subplot(2, 4, sess); hold on
    for perm = 1 : 3
        plot(out_samp(perm).tt, lev_ctrl(count).choice.sig.rho(:, perm, sess), 'k')
        plot(out_samp(perm).tt, lev_ctrl(count).choice.sig.rho_part(:, perm, sess), 'r')
        plot(out_samp(perm).tt, lev_ctrl(count).choice.sig.rho_part_other(:, perm, sess), 'g')
    end
end

figure; hold on
for sess = 1 : 8
    subplot(2, 4, sess); hold on
    for perm = 1 : 3
        plot(out_samp(perm).tt, lev_ctrl(count).choice.perm.rho(:, perm, sess), 'k')
        plot(out_samp(perm).tt, lev_ctrl(count).choice.perm.rho_part(:, perm, sess), 'r')
        plot(out_samp(perm).tt, lev_ctrl(count).choice.perm.rho_part_other(:, perm, sess), 'g')
    end
end



figure; hold on
for perm = 1 : 10
plot(out_samp(perm).tt, out_samp(perm).rho_RT, 'k')
plot(out_samp(perm).tt, out_samp(perm).rho_RT_partial, 'r')
end

figure; hold on
for perm = 1 : 10
plot(out_samp(perm).tt, out_samp(perm).rho_choice, 'k')
plot(out_samp(perm).tt, out_samp(perm).rho_choice_partial, 'r')
end

out_perm


figure; hold on
for perm = 1 : 10
plot(out_samp(perm).tt, out_perm(perm).rho_RT, 'k')
plot(out_samp(perm).tt, out_perm(perm).rho_RT_partial, 'r')
end

figure; hold on
for perm = 1 : 10
plot(out_samp(perm).tt, out_perm(perm).rho_choice, 'k')
plot(out_samp(perm).tt, out_perm(perm).rho_choice_partial, 'r')
end










