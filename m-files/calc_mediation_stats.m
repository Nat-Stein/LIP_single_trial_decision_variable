



medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));

mediatedDimsA = [52 19 27]; dimTin = 27; % ramp, PC1, Tin



%% =======================================================================
%% Compare observed value to distribution of values in shuffled trial assignments
proj_type = 'PC1weights'; % 'PC1weights', 'uniform', 'randn', 'gauss'
load(fullfile(saveLoc, ['mediation_randProj_' proj_type]))

inclSh = 1 : 1000;
med_stats = [];
for mDim = mediatedDimsA
    
    m = medA.leverage(mDim).choice.unmediated.m;
    mRT = medA.leverage(mDim).RT.unmediated.m;
    
    for tp = 1 : length(m)
        
        % Choice
        dat = m(tp);
        null_distrib = squeeze(nanmean(med_randProj.rho_choice(:, inclSh, tp), 1));
        
        med_stats.choice.pDirect{mDim}(tp) = sum(null_distrib > dat) / length(inclSh);
        
        % RT
        dat = mRT(tp);
        null_distrib = squeeze(nanmean(med_randProj.rho_RT(:, inclSh, tp), 1));
        med_stats.RT.pDirect{mDim}(tp) = sum(null_distrib < dat) / length(inclSh);
    end 
end

figure; hold on
subplot(1, 2, 1); hold on
for mDim = mediatedDimsA
    plot(med_shuffle.t, med_stats.choice.pDirect{mDim})
end
% p < 0.05 starting at 270ms for all
% starting 230ms for ramp, 190ms for PC1 and 270ms for TinCon
subplot(1, 2, 2); hold on
for mDim = mediatedDimsA
    plot(med_shuffle.t, med_stats.RT.pDirect{mDim})
end
% p < 0.05 starting at 280ms
% starting 310ms for ramp, 300ms for PC1 and 180ms for TinCon
save(fullfile(saveLoc, 'med_stats_PC1shuffle'), 'med_stats')



%% =======================================================================
%% Compare observed value to distribution of values in shuffled trial assignments
load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se_shuffle'));
inclSh = 1 : 500;

med_stats = [];
for mDim = mediatedDimsA
    
    m = medA.leverage(mDim).choice.unmediated.m;
    mRT = medA.leverage(mDim).RT.unmediated.m;
    
    for tp = 1 : length(m)
        
        % Choice
        dat = m(tp);
        null_distrib = squeeze(med_shuffle.rho_choice(inclSh, mDim, tp));
        
        med_stats.choice.pDirect{mDim}(tp) = sum(null_distrib > dat) / length(inclSh);
        
        % RT
        dat = mRT(tp);
        null_distrib = squeeze(med_shuffle.rho_RT(inclSh, mDim, tp));
        med_stats.RT.pDirect{mDim}(tp) = sum(null_distrib < dat) / length(inclSh);
    end 
end

figure; hold on
subplot(1, 2, 1); hold on
for mDim = mediatedDimsA
    plot(med_shuffle.t, med_stats.choice.pDirect{mDim})
end
% p < 0.05 starting at 240ms for all
% starting 190ms for ramp, 170ms for PC1 and 240ms for TinCon
subplot(1, 2, 2); hold on
for mDim = mediatedDimsA
    plot(med_shuffle.t, med_stats.RT.pDirect{mDim})
end
% p < 0.05 starting at 280ms
% starting 280ms for ramp, 260ms for PC1 and 180ms for TinCon
save(fullfile(saveLoc, 'med_stats'), 'med_stats')







