% Correlation with behavior and mediation of the effects


%% Normalize activity in each session based on minimum and maximum amplitude during decision epoch (same as in singleTrial_....m)

normDim_sl = [];
normDim_rl = [];

dims = [1 2 7 20 31 32 63]; %dims = 20;
fw = ones(1, 50)/50;
tsl = find(par.tsl >= 0 & par.tsl <=0.6); 
trl = find(par.trl >= -0.6 & par.trl <=0); 
for dim = dims
    spkSL = [];
    spkRL = [];
    for sess = 1 : 8
        par = par_sess{sess};
        sc = 1; sn = 1; % if dim == 20; sc = pcTin(sess); sn = pcTin_sign(sess); end
        
        newSL = sn * squeeze(allDim_sl{dim}{sc, sess});
        newRL = sn * squeeze(allDim_rl{dim}{sc, sess});
        
        mm = [];
        for ch = 0 : 1
            c0 = find(par.dataMat(:, par.idx.cho_trg) == ch);
            add = conv(squeeze(nanmean(newSL(c0, :), 1)), fw, 'same');
            mm = [mm add(1, tsl)];
            add = conv(squeeze(nanmean(newRL(c0, :), 1)), fw, 'same');
            mm = [mm add(1, trl)];
        end
        maxM = max(mm);
        minM = min(mm);
        exc = maxM - minM;
        
        normDim_paraS{dim}{sc, sess}.mini = minM;
        normDim_paraS{dim}{sc, sess}.maxi = maxM;
        normDim_paraS{dim}{sc, sess}.excur = exc;
        normDim_paraS{dim}{sc, sess}.mm = mm;
        
        normSL = (newSL - minM) ./ exc;
        normRL = (newRL - minM) ./ exc;
        
        normDim_sl{dim}{sc, sess} = normSL;
        normDim_rl{dim}{sc, sess} = normRL;
        
    end
end

% save('E:\Matalb analyses\norm_dim_sl_rl', 'normDim_sl', 'normDim_rl', 'normDim_para', 'par_sess', '-v7.3')

%% amp vs choice & RT normalized by max activity on Tin trials - STILL IN USE - YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
% Prepare data
he=1;
bWidth = 0.08;
whichTs = -0.1 : 0.005 : 0.55; whichTs = round(whichTs, 3);
dataMat = par_sess{sess}.dataMat;
idx = par_sess{sess}.idx;
cohs = unique(dataMat(:, idx.sig_coh)); cohs = cohs(find(abs(cohs) <= 0.2))';
miniRT = 0.67;
maxiRT = 2;
t_CBLate = [.525 .575]; 
tComp = find(par.tsl > t_CBLate(1) & par.tsl < t_CBLate(2));

useDims = [1 2 7 20 31 32 63];

% This has nothing to do with residuals, right? YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
clear resDimT_norm resDim_norm
for sess = 1 : 8
    dataMat = par_sess{sess}.dataMat;
    rt = [];
    choice = [];
    sigcoh = [];
    for dim = useDims
        sc = 1; sn = 1; % if dim == 20; sc = pcTin(sess); sn = pcTin_sign(sess); end
        compAmp = [];
        
        for tt = 1 : length(whichTs)
            amp = [];
            trls = find(dataMat(:, idx.task) == he &...
                dataMat(:, idx.effector) == he &...
                dataMat(:, idx.rt) > miniRT*1000 &...
                dataMat(:, idx.rt) < maxiRT*1000);
            
            tpts = find(par.tsl >= whichTs(tt)-bWidth/2 & par.tsl <= whichTs(tt)+bWidth/2);
            dat = nanmean(normDim_sl{dim}{sc, sess}(trls, tpts), 2);
            amp = [amp dat'];
            if tt == 1
                dat = nanmean(normDim_sl{dim}{sc, sess}(trls, tComp), 2);
                compAmp = [compAmp dat'];
                if dim == useDims(1)
                    RT = dataMat(trls, idx.rt);
                    rt = [rt RT'];
                    c = 1 - dataMat(trls, idx.cho_trg);
                    choice = [choice c'];
                    c = dataMat(trls, idx.sig_coh);
                    sigcoh = [sigcoh c'];
                end
            end
            resDimT_norm.amp{dim, sess}(:, tt) = amp;
            
        end
        resDimT_norm.compAmp{dim, sess} = compAmp;
    end
    resDim_norm.rt{sess} = rt;
    resDim_norm.choice{sess} = choice;
    resDim_norm.sigcoh{sess} = sigcoh;
end

% Choice ----------------------------------------------------------------
% clear dim_predChoPool dim_predXPool
maxCoh = 0.2;
medDims = useDims;
for dim = useDims
    for dimM = medDims
        disp([num2str(dim) ' vs ' num2str(dimM)])
        
        for tt = 1 : length(whichTs)
            amp = [];
            ampM = [];
            allCoh = [];
            choice = [];
            % session = [];
            for sess = 1 : 8
                
                trls = find(resDim_norm.sigcoh{sess} < maxCoh);
                
                addAmp = resDimT_norm.amp{dim, sess}(trls, tt)';
                amp = [amp addAmp];
                addAmp = resDimT_norm.compAmp{dimM, sess}(1, trls);
                ampM = [ampM addAmp];
                c = resDim_norm.sigcoh{sess}(1, trls);
                allCoh = [allCoh c];
                c = resDim_norm.choice{sess}(1, trls);
                choice = [choice c];
                % session = [session sess * ones(1, length(c))];
            end
            cohs = unique(allCoh);
            cohMat = [];
            for c = 1 : length(cohs)
                cohMat(c, :) = (allCoh == cohs(c));
            end
            %             subjMat = [];
            %             for c = unique(session)
            %                 subjMat(c, :) = (session == c);
            %             end
            
            if dim == dimM
                
                % Choice - SL: raw, coh
                utrls = ~isnan(amp) & ~isnan(choice);
                y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
                % x = [amp(utrls)' allCoh(utrls)' subjMat(:, utrls)'];
                x = [amp(utrls)'];
                
                [b, dev, stat] = glmfit(x,y,...
                    'binomial', 'link', 'logit');
                
                dim_predChoPool.ChoRes_normB_d1_sess{dim}(tt, :) = stat.beta;       % new
                dim_predChoPool.ChoRes_normp_d1_sess{dim}(tt, :) = stat.p;
                dim_predChoPool.ChoRes_normT_d1_sess{dim}(tt, :) = stat.t;
                dim_predChoPool.ChoRes_normSE_d1_sess{dim}(tt, :) = stat.se;        % new
                dim_predChoPool.ChoRes_normVarEx_d1_sess{dim}(tt, :) = 1-var(stat.resid)./var(y);
                
                % Choice - SL: raw, coh
                utrls = ~isnan(amp) & ~isnan(choice) & ~isnan(allCoh);
                y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
                % x = [amp(utrls)' allCoh(utrls)' subjMat(:, utrls)'];
                x = [amp(utrls)' allCoh(utrls)'];
                
                [b, dev, stat] = glmfit(x,y,...
                    'binomial', 'link', 'logit');
                
                dim_predChoPool.ChoRes_normB_coh_d1_sess{dim}(tt, :) = stat.beta;       % WE USE THIS ONE
                dim_predChoPool.ChoRes_normp_coh_d1_sess{dim}(tt, :) = stat.p;
                dim_predChoPool.ChoRes_normT_coh_d1_sess{dim}(tt, :) = stat.t;
                dim_predChoPool.ChoRes_normSE_coh_d1_sess{dim}(tt, :) = stat.se;        % WE USE THIS ONE
                dim_predChoPool.ChoRes_normVarEx_coh_d1_sess{dim}(tt, :) = 1-var(stat.resid)./var(y);
                
                
                
                % Choice - SL: raw. cohMat
                utrls = ~isnan(amp) & ~isnan(choice) & ~isnan(allCoh);
                y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
                % x = [amp(utrls)' cohMat(:, utrls)' subjMat(:, utrls)'];
                x = [amp(utrls)' cohMat(:, utrls)'];
                
                [b, dev, stat] = glmfit(x,y,...
                    'binomial', 'link', 'logit');
                
                dim_predChoPool.ChoRes_normB_cohMat_d1_sess{dim}(tt, :) = stat.beta;
                dim_predChoPool.ChoRes_normp_cohMat_d1_sess{dim}(tt, :) = stat.p;
                dim_predChoPool.ChoRes_normT_cohMat_d1_sess{dim}(tt, :) = stat.t;
                dim_predChoPool.ChoRes_normSE_cohMat_d1_sess{dim}(tt, :) = stat.se;
                dim_predChoPool.ChoRes_normVarEx_cohMat_d1_sess{dim}(tt, :) = 1-var(stat.resid)./var(y);
                
            end
            
            % Choice - SL: partial, coh
            utrls = ~isnan(amp) & ~isnan(ampM) & ~isnan(choice);
            y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
            % x = [amp(utrls)' allCoh(utrls)' ampM(utrls)' subjMat(:, utrls)'];
            x = [amp(utrls)' ampM(utrls)'];
            
            [b, dev, stat] = glmfit(x,y,...
                'binomial', 'link', 'logit');
            
            dim_predChoPool.Cho_xMedRes_norm_d1_B_sess{dim, dimM}(tt, :) = stat.beta;        % WE USE THIS ONE
            dim_predChoPool.Cho_xMedRes_norm_d1_p_sess{dim, dimM}(tt, :) = stat.p;
            dim_predChoPool.Cho_xMedRes_norm_d1_T_sess{dim, dimM}(tt, :) = stat.t;
            dim_predChoPool.Cho_xMedRes_norm_d1_SE_sess{dim, dimM}(tt, :) = stat.se;         % WE USE THIS ONE
            dim_predChoPool.Cho_xMedRes_norm_d1_VarEx_sess{dim, dimM}(tt, :) = 1-var(stat.resid)./var(y);
            
            
            % Choice - SL: partial, coh
            utrls = ~isnan(amp) & ~isnan(ampM) & ~isnan(choice) & ~isnan(allCoh);
            y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
            % x = [amp(utrls)' allCoh(utrls)' ampM(utrls)' subjMat(:, utrls)'];
            x = [amp(utrls)' allCoh(utrls)' ampM(utrls)'];
            
            [b, dev, stat] = glmfit(x,y,...
                'binomial', 'link', 'logit');
            
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess{dim, dimM}(tt, :) = stat.beta;        % WE USE THIS ONE
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_p_sess{dim, dimM}(tt, :) = stat.p;
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_T_sess{dim, dimM}(tt, :) = stat.t;
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_SE_sess{dim, dimM}(tt, :) = stat.se;         % WE USE THIS ONE
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_VarEx_sess{dim, dimM}(tt, :) = 1-var(stat.resid)./var(y);
            
            % Choice - SL: partial, cohMat
            utrls = ~isnan(amp) & ~isnan(ampM) & ~isnan(choice) & ~isnan(allCoh);
            y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
            % x = [amp(utrls)' cohMat(:, utrls)' ampM(utrls)' subjMat(:, utrls)'];
            x = [amp(utrls)' cohMat(:, utrls)' ampM(utrls)'];
            
            [b, dev, stat] = glmfit(x,y,...
                'binomial', 'link', 'logit');
            
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_B_sess{dim, dimM}(tt, :) = stat.beta;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_p_sess{dim, dimM}(tt, :) = stat.p;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_T_sess{dim, dimM}(tt, :) = stat.t;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_SE_sess{dim, dimM}(tt, :) = stat.se;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_VarEx_sess{dim, dimM}(tt, :) = 1-var(stat.resid)./var(y);
        end
    end
end


% RT ----------------------------------------------------------------
medDims = useDims;
for cho = 1 : 2
    for dim = useDims
        for dimM = medDims
            disp([num2str(dim) ' vs ' num2str(dimM)])
            
            for tt = 1 : length(whichTs)
                amp = [];
                ampM = [];
                allCoh = [];
                rt = [];
                % session = [];
                for sess = 1 : 8
                    trls = find(resDim_norm.choice{sess} == 2-cho &...
                        resDim_norm.sigcoh{sess} < maxCoh);
                    addAmp = resDimT_norm.amp{dim, sess}(trls, tt)';
                    amp = [amp addAmp];
                    addAmp = resDimT_norm.compAmp{dimM, sess}(1, trls);
                    ampM = [ampM addAmp];
                    c = resDim_norm.sigcoh{sess}(1, trls);
                    allCoh = [allCoh c];
                    c = resDim_norm.rt{sess}(1, trls);
                    rt = [rt c];
                    % session = [session sess * ones(1, length(c))];
                end
                cohs = unique(allCoh);
                cohMat = [];
                for c = 1 : length(cohs)
                    cohMat(c, :) = (allCoh == cohs(c));
                end
                %             subjMat = [];
                %             for c = unique(session)
                %                 subjMat(c, :) = (session == c);
                %             end
                
                if dim == dimM
                    % RT - SL normalized (not zscored)
                    utrls = ~isnan(amp) & ~isnan(rt);
                    y = rt(utrls)';
                    x = amp(utrls)';
                    
                    [R,P] = corrcoef(x, y, 'rows', 'pairwise');
                    
                    dim_predXPool.pre_normRTr{dim, cho}(tt) = R(1, 2);
                    dim_predXPool.pre_normRTp{dim, cho}(tt) = P(1, 2);
                    dim_predXPool.pre_normRT_N(cho) = length(x);            % WE USE THIS ONE
                    
                    % RT - SL normalized with coherence indicator (not zscored)
                    utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt);
                    y = rt(utrls)';
                    x = amp(utrls)';
                    z = allCoh(utrls)';
                    
                    [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                    
                    dim_predXPool.pre_normCohRTr{dim, cho}(tt) = rho;       % WE USE THIS ONE
                    dim_predXPool.pre_normCohRTp{dim, cho}(tt) = p;
                    
                    % RT - SL normalized with coherence indicator (not zscored)
                    utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt);
                    y = rt(utrls)';
                    x = amp(utrls)';
                    z = cohMat(:, utrls)';
                    
                    [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                    
                    dim_predXPool.pre_normCohMatRTr{dim, cho}(tt) = rho;
                    dim_predXPool.pre_normCohMatRTp{dim, cho}(tt) = p;
                    
                    
                end
                % RT - SL normalized (not zscored) - mediated
                utrls = ~isnan(amp) & ~isnan(rt) & ~isnan(ampM);
                y = rt(utrls)';
                x = amp(utrls)';
                z = ampM(utrls)';
                
                
                [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                
                dim_predXPool.preX_normRTr{dim, dimM, cho}(tt) = rho;
                dim_predXPool.preX_normRTp{dim, dimM, cho}(tt) = p;
                
                % RT - SL normalized with coherence indicator (not zscored) - mediated
                utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt) & ~isnan(ampM);
                y = rt(utrls)';
                x = amp(utrls)';
                z = [allCoh(utrls)' ampM(utrls)'];
                
                [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                
                dim_predXPool.preX_normCohRTr{dim, dimM, cho}(tt) = rho;    % WE USE THIS ONE
                dim_predXPool.preX_normCohRTp{dim, dimM, cho}(tt) = p;
                
                % RT - SL normalized with coherence indicator (not zscored) - mediated
                utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt) & ~isnan(ampM);
                y = rt(utrls)';
                x = amp(utrls)';
                z = [cohMat(:, utrls)' ampM(utrls)'];
                
                [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                
                dim_predXPool.preX_normCohMatRTr{dim, dimM, cho}(tt) = rho;
                dim_predXPool.preX_normCohMatRTp{dim, dimM, cho}(tt) = p;
                
            end
        end
    end
end

% save('E:\Matalb analyses\dim_pred_pool_230329', 'dim_predChoPool', 'dim_predXPool', 'par_sess', '-v7.3')
% save('E:\Matalb analyses\dim_pred_pool_230329_080', 'dim_predChoPool', 'dim_predXPool', 'par_sess', '-v7.3')

%% Plot Amp-vs-RT for data pooled across sessions - neuron classes (+med of PC1/ramp) with standard error - THIS ONE - THIS ONE - THIS ONE - THIS ONE - THIS ONE - THIS ONE - THIS ONE - THIS ONE 

% load('E:\Matalb analyses\dim_pred_pool_230329_080')
% load('E:\Matalb analyses\dim_pred_pool_230329')

xLims = [0.2 0.56]; 
clear yLimsRT yLimsCho
dim = 20; yLimsRT{dim} = [-0.3 0.15]; yLimsCho{dim} = [-1 2];
dim = 63; yLimsRT{dim} = [-0.3 0.15]; yLimsCho{dim} = [-1 3];

compTP = find(round(whichTs, 3) == 0.55);
plotSess = [1 : 8];
barCols = {[233 149 53]/255, [76 0 153]/255, [0 153 76]/255, [233 53 149]/255,...
    [0 149 233]/255, [9 237 104]/255, [76 0 153]/255, [149 53 149]/255, [0 153 76]/255};

pltDims = [1 2 31 32]; mediatedDims = [20 63];

pltDim = dimName; pltDim{1} = 'Tin-c'; pltDim{2} = 'Tin-i'; pltDim{20} = 'PC1'; pltDim{63} = 'rampD'; pltDim{69} = 'Tin w_r_a_m_p';
rows = 2; colms = 2;
plotTraw = find(whichTs >=0.2 & whichTs <=0.45);

% RT correlations based on data normalized across sessions using the max
% amplitude and excursion in decision period
% rawB = dim_predXPool.pre_normRTr; 
% medR = dim_predXPool.preX_normRTr; 
rawB = dim_predXPool.pre_normCohRTr; % Which one to use? This one considers coherence and should be it, right?
medR = dim_predXPool.preX_normCohRTr;
% N = dim_predXPool.pre_normRT_N; 
% medSE = dim_predXPool.pre_normCohRTr;
% rawSE = sqrt((1 - rawB.^2)./(N(cho)-2));
% -----------------------
% Leverage on choice
rawB_cho = dim_predChoPool.ChoRes_normB_d1_sess;
rawSE_cho = dim_predChoPool.ChoRes_normSE_d1_sess;
medR_cho = dim_predChoPool.Cho_xMedRes_norm_d1_B_sess;
medSE_cho = dim_predChoPool.Cho_xMedRes_norm_d1_SE_sess;
% rawB_cho = dim_predChoPool.ChoRes_normB_coh_d1_sess;
% rawSE_cho = dim_predChoPool.ChoRes_normSE_coh_d1_sess;
% medR_cho = dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess;
% medSE_cho = dim_predChoPool.Cho_xMedRes_norm_coh_d1_SE_sess;

% rawB_cho = dim_predChoPool.ChoRes_normB_cohMat_d1_sess;
% rawSE_cho = dim_predChoPool.ChoRes_normSE_cohMat_d1_sess;
% medR_cho = dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_B_sess;
% medSE_cho = dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_SE_sess;

saveFig = 0;

sFont = 6;
compTP = find(whichTs == 0.55);
for mDim = mediatedDims
    
    plotDims = [pltDims mDim];
    figure('Position', [200 10 700 700]); hold on
    
    % Choice =============================================================
    % Leverage on choice: neuron classes
    subplot(rows, colms, 1); hold on
    for d1 = 1 : length(plotDims)
        
        dim = plotDims(d1);
        pltMeanRaw = squeeze(rawB_cho{dim}(:, 2));
        pltSERaw = squeeze(rawSE_cho{dim}(:, 2));
        
        % All dimensions raw
        plotCol = barCols{d1};
        s = shadedErrorBar(whichTs(plotTraw), pltMeanRaw(plotTraw), pltSERaw(plotTraw));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = 3;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        pltMeanRaw = squeeze(rawB_cho{dim}(compTP, 2));
        pltSERaw = squeeze(rawSE_cho{dim}(compTP, 2));
        plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
        line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    end
    title('Raw')
    ylabel('Leverage on choice (slope beta)')
    xlabel('Time after motion onset [s]')
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimsCho{mDim})
    % ylim(max(abs(yLims))*[-1 1]);
    set(gca,'TickDir','out', 'FontSize', sFont);
    
    % PC1 and Ramp dimension mediated by all neuron classes
    dim = mDim;
    subplot(rows, colms, 3); hold on; title(pltDim{dim})
    pltMeanRaw = squeeze(rawB_cho{dim}(:, 2));
    pltSERaw = squeeze(rawSE_cho{dim}(:, 2));
    plotCol = 'k';
    s = shadedErrorBar(whichTs(plotTraw), pltMeanRaw(plotTraw), pltSERaw(plotTraw));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = 3;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    pltMeanRaw = squeeze(rawB_cho{dimM}(compTP, 2));
    pltSERaw = squeeze(rawSE_cho{dimM}(compTP, 2));
    plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
    line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    
    for d1 = 1 : length(plotDims)
        dimM = plotDims(d1);
        plotCol = barCols{d1};
        pltMean = squeeze(medR_cho{dim, dimM}(:, 2));
        pltSE = squeeze(medSE_cho{dim, dimM}(:, 2));
        s = shadedErrorBar(whichTs(plotTraw), pltMean(plotTraw), pltSE(plotTraw));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = 2;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        pltMeanRaw = squeeze(rawB_cho{dimM}(compTP, 2));
        pltSERaw = squeeze(rawSE_cho{dimM}(compTP, 2));
        plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
        line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    end
    ylabel('Leverage on choice (slope beta)')
    xlabel('Time from motion onset [s]')
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimsCho{mDim}); if mDim == 2 || mDim == 32; ylim(-flip(yLimsCho)); end
    set(gca,'TickDir','out', 'FontSize', sFont);
    % legend('Raw', 'Tin-c', 'Tin-i', 'Min-c', 'Min-i', 'PC1', 'ramp')
    
    
    % RT =============================================================
    % RT correlations: neuron classes
    subplot(rows, colms, 2); hold on
    for d1 = 1 : length(plotDims)
        
        dim = plotDims(d1);
        pltMeanRaw = squeeze(rawB{dim, cho});
        pltSERaw = sqrt((1 - pltMeanRaw.^2)./(N(cho)-2));
        
        % All dimensions raw
        plotCol = barCols{d1};
        s = shadedErrorBar(whichTs(plotTraw), pltMeanRaw(plotTraw), pltSERaw(plotTraw));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = 3;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        pltMeanRaw = squeeze(rawB{dim, cho}(compTP));
        pltSERaw = sqrt((1 - pltMeanRaw.^2)./(N(cho)-2));
        plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
        line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    end
    title('Raw')
    ylabel('Correltaion with RT (r)')
    xlabel('Time after motion onset [s]')
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimsRT{mDim})
    % ylim(max(abs(yLims))*[-1 1]);
    set(gca,'TickDir','out', 'FontSize', sFont);
    
    % PC1 and Ramp dimension mediated by all neuron classes
    legs = []; l = 0; 
    dim = mDim;
    subplot(rows, colms, 4); hold on; title(pltDim{dim})
    pltMeanRaw = squeeze(rawB{dim, cho});
    pltSERaw = sqrt((1 - pltMeanRaw.^2)./(N(cho)-2));
    plotCol = 'k';
    s = shadedErrorBar(whichTs(plotTraw), pltMeanRaw(plotTraw), pltSERaw(plotTraw));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = 3;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    l = l + 1; legs{l} = 'Raw';
    
    for d1 = 1 : length(plotDims)
        dimM = plotDims(d1);
        plotCol = barCols{d1};
        pltMean = squeeze(medR{dim, dimM, cho});
        pltSE = sqrt((1 - pltMean.^2)./(N(cho)-2));
        % pltSE = squeeze(medSE{dim, dimM}(:, 2));
        s = shadedErrorBar(whichTs(plotTraw), pltMean(plotTraw), pltSE(plotTraw));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = 2;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        l = l + 1; legs{l} = dimName{dimM};
    end
    for d1 = 1 : length(plotDims)
        dimM = plotDims(d1);
        plotCol = barCols{d1};
        pltMean = squeeze(medR{dim, dimM, cho});
        pltSE = sqrt((1 - pltMean.^2)./(N(cho)-2));
        
        pltMeanRaw = squeeze(rawB{dimM, cho}(compTP));
        pltSERaw = sqrt((1 - pltMeanRaw.^2)./(N(cho)-2));
        plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
        line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    end
    ylabel('Correltaion with RT (r)')
    xlabel('Time from motion onset [s]')
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimsRT{mDim}); if mDim == 2 || mDim == 32; ylim(-flip(yLimsRT)); end
    set(gca,'TickDir','out', 'FontSize', sFont);
    legend(legs)
    
    if saveFig == 1
        ff = fullfile('C:\Users\shadlenlab\Dropbox\LIP single trial\Figures\behavior correlation',...
            [dimName{dimM} '_RT_Cho_med_classes_allSess']);
        saveas(gca, [ff '.fig'])
        saveas(gca, [ff '.pdf'])
        saveas(gca, [ff '.png'])
    end
    
end


%% Remake other paper figure


pltDims = [20 63 1]; mediatedDims = [20 63 1];
%rows = length(mediatedDims); colms = 2; 
rows = 2; colms = length(mediatedDims); 

saveFig = 0;

sFont = 6;
compTP = find(whichTs == 0.55);
% figure('Position', [200 10 700 1000]); hold on
figure('Position', [200 10 1000 700]); hold on
ro = 0; 
for mDim = mediatedDims
    
    ro = ro + 1;
    plotDims = pltDims;
    % Choice =============================================================
    % Leverage on choice: neuron classes
    % subplot(rows, colms, (ro-1)*colms + 1); hold on; title(pltDim{mDim})
    subplot(rows, colms, ro); hold on; title(pltDim{mDim})
    dim = mDim;
    pltMeanRaw = squeeze(rawB_cho{dim}(:, 2));
    pltSERaw = squeeze(rawSE_cho{dim}(:, 2));
    plotCol = 'k';
    s = shadedErrorBar(whichTs(plotTraw), pltMeanRaw(plotTraw), pltSERaw(plotTraw));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = 3;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    pltMeanRaw = squeeze(rawB_cho{dimM}(compTP, 2));
    pltSERaw = squeeze(rawSE_cho{dimM}(compTP, 2));
    plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
    line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    
    for d1 = 1 : length(plotDims)
        dimM = plotDims(d1);
        plotCol = barCols{d1};
        pltMean = squeeze(medR_cho{dim, dimM}(:, 2));
        pltSE = squeeze(medSE_cho{dim, dimM}(:, 2));
        s = shadedErrorBar(whichTs(plotTraw), pltMean(plotTraw), pltSE(plotTraw));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = 2;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        pltMeanRaw = squeeze(rawB_cho{dimM}(compTP, 2));
        pltSERaw = squeeze(rawSE_cho{dimM}(compTP, 2));
        plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
        line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    end
    ylabel('Leverage on choice (slope beta)')
    xlabel('Time from motion onset [s]')
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    % ylim(yLimsCho{mDim}); if mDim == 2 || mDim == 32; ylim(-flip(yLimsCho)); end
    set(gca,'TickDir','out', 'FontSize', sFont);
    % legend('Raw', 'Tin-c', 'Tin-i', 'Min-c', 'Min-i', 'PC1', 'ramp')
    
    
    % RT =============================================================
    % RT correlations: neuron classes
    
    % PC1 and Ramp dimension mediated by all neuron classes
    legs = []; l = 0; 
    dim = mDim;
    % subplot(rows, colms, (ro-1)*colms + 2); hold on; title(pltDim{mDim})
    subplot(rows, colms, colms + ro); hold on; title(pltDim{mDim})
    pltMeanRaw = squeeze(rawB{dim, cho});
    pltSERaw = sqrt((1 - pltMeanRaw.^2)./(N(cho)-2));
    plotCol = 'k';
    s = shadedErrorBar(whichTs(plotTraw), pltMeanRaw(plotTraw), pltSERaw(plotTraw));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = 3;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    l = l + 1; legs{l} = 'Raw';
    
    for d1 = 1 : length(plotDims)
        dimM = plotDims(d1);
        plotCol = barCols{d1};
        pltMean = squeeze(medR{dim, dimM, cho});
        pltSE = sqrt((1 - pltMean.^2)./(N(cho)-2));
        % pltSE = squeeze(medSE{dim, dimM}(:, 2));
        s = shadedErrorBar(whichTs(plotTraw), pltMean(plotTraw), pltSE(plotTraw));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = 2;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        l = l + 1; legs{l} = dimName{dimM};
    end
    pltMeanRaw = squeeze(rawB{dim, cho}(compTP));
    pltSERaw = sqrt((1 - pltMeanRaw.^2)./(N(cho)-2));
    plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
    line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    for d1 = 1 : length(plotDims)
        dimM = plotDims(d1);
        plotCol = barCols{d1};
        pltMean = squeeze(medR{dim, dimM, cho});
        pltSE = sqrt((1 - pltMean.^2)./(N(cho)-2));
        
        pltMeanRaw = squeeze(rawB{dimM, cho}(compTP));
        pltSERaw = sqrt((1 - pltMeanRaw.^2)./(N(cho)-2));
        plot(whichTs(compTP), pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
        line(whichTs(compTP)*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    end
    ylabel('Correltaion with RT (r)')
    xlabel('Time from motion onset [s]')
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    %ylim(yLimsRT{mDim}); if mDim == 2 || mDim == 32; ylim(-flip(yLimsRT)); end
    set(gca,'TickDir','out', 'FontSize', sFont);
    legend(legs)
    
    %     if saveFig == 1
    %         ff = fullfile('C:\Users\shadlenlab\Dropbox\LIP single trial\Figures\behavior correlation',...
    %             [dimName{dimM} '_RT_Cho_med_classes_allSess']);
    %         saveas(gca, [ff '.fig'])
    %         saveas(gca, [ff '.pdf'])
    %         saveas(gca, [ff '.png'])
    %     end
    
end





%% OLD?
%% Compute residuals within each motion strength
% Prepare data
% bWidth = 0.05;
% whichTs = -0.1 : 0.005 : 0.55; whichTs = round(whichTs, 3);
% dataMat = par_sess{sess}.dataMat;
% idx = par_sess{sess}.idx;
% cohs = unique(dataMat(:, idx.sig_coh)); cohs = cohs(find(abs(cohs) <= 0.2))';
for sess = 1 : 8
    dataMat = par_sess{sess}.dataMat;
    rt = [];
    choice = [];
    sigcoh = [];
    for dim = useDims
        sc = 1; sn = 1; % if dim == 20; sc = pcTin(sess); sn = pcTin_sign(sess); end
        compAmp = [];
        
        for tt = 1 : length(whichTs)
            amp = [];
            trls = find(dataMat(:, idx.task) == he &...
                dataMat(:, idx.effector) == he &...
                dataMat(:, idx.rt) > miniRT &...
                dataMat(:, idx.rt) < maxiRT);
            
            tpts = find(par.tsl >= whichTs(tt)-bWidth/2 & par.tsl <= whichTs(tt)+bWidth/2);
            dat = nanmean(normDim_sl{dim}{sc, sess}(trls, tpts), 2);
            amp = [amp dat'];
            if tt == 1
                dat = nanmean(normDim_sl{dim}{sc, sess}(trls, tComp), 2);
                compAmp = [compAmp dat'];
                if dim == useDims(1)
                    RT = dataMat(trls, idx.rt);
                    rt = [rt RT'];
                    c = 1 - dataMat(trls, idx.cho_trg);
                    choice = [choice c'];
                    c = dataMat(trls, idx.sig_coh);
                    sigcoh = [sigcoh c'];
                end
            end
            resDimT_norm{dim, sess}(:, tt) = amp;
            
        end
        resDimT_norm_compAmp{dim, sess} = compAmp;
    end
    resDim_norm_rt{sess} = rt;
    resDim_norm_choice{sess} = choice;
    resDim_norm_sigcoh{sess} = sigcoh;
end

% Choice ----------------------------------------------------------------
medDims = useDims;
for dim = useDims
    for dimM = medDims
        disp([num2str(dim) ' vs ' num2str(dimM)])
        
        for tt = 1 : length(whichTs)
            amp = [];
            ampM = [];
            allCoh = [];
            choice = [];
            % session = [];
            for sess = 1 : 8
                
                addAmp = resDimT_norm{dim, sess}(:, tt)';
                amp = [amp addAmp];
                addAmp = resDimT_norm_compAmp{dimM, sess};
                ampM = [ampM addAmp];
                c = resDim_norm_sigcoh{sess};
                allCoh = [allCoh c];
                c = resDim_norm_choice{sess};
                choice = [choice c];
                % session = [session sess * ones(1, length(c))];
            end
            cohs = unique(allCoh);
            cohMat = [];
            for c = 1 : length(cohs)
                cohMat(c, :) = (allCoh == cohs(c));
            end
            %             subjMat = [];
            %             for c = unique(session)
            %                 subjMat(c, :) = (session == c);
            %             end
            
            if dim == dimM
                % Choice - zscored SL: raw, coh
                utrls = ~isnan(amp) & ~isnan(choice) & ~isnan(allCoh);
                y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
                % x = [amp(utrls)' allCoh(utrls)' subjMat(:, utrls)'];
                x = [amp(utrls)' allCoh(utrls)'];
                
                [b, dev, stat] = glmfit(x,y,...
                    'binomial', 'link', 'logit');
                
                dim_predChoPool.ChoRes_normB_coh_d1_sess{dim}(tt, :) = stat.beta;
                dim_predChoPool.ChoRes_normp_coh_d1_sess{dim}(tt, :) = stat.p;
                dim_predChoPool.ChoRes_normT_coh_d1_sess{dim}(tt, :) = stat.t;
                dim_predChoPool.ChoRes_normSE_coh_d1_sess{dim}(tt, :) = stat.se;
                dim_predChoPool.ChoRes_normVarEx_coh_d1_sess{dim}(tt, :) = 1-var(stat.resid)./var(y);
                
                
                
                % Choice - zscored SL: raw. cohMat
                utrls = ~isnan(amp) & ~isnan(choice) & ~isnan(allCoh);
                y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
                % x = [amp(utrls)' cohMat(:, utrls)' subjMat(:, utrls)'];
                x = [amp(utrls)' cohMat(:, utrls)'];
                
                [b, dev, stat] = glmfit(x,y,...
                    'binomial', 'link', 'logit');
                
                dim_predChoPool.ChoRes_normB_cohMat_d1_sess{dim}(tt, :) = stat.beta;
                dim_predChoPool.ChoRes_normp_cohMat_d1_sess{dim}(tt, :) = stat.p;
                dim_predChoPool.ChoRes_normT_cohMat_d1_sess{dim}(tt, :) = stat.t;
                dim_predChoPool.ChoRes_normSE_cohMat_d1_sess{dim}(tt, :) = stat.se;
                dim_predChoPool.ChoRes_normVarEx_cohMat_d1_sess{dim}(tt, :) = 1-var(stat.resid)./var(y);
                
            end
            
            % Choice - zscored SL: partial, coh
            utrls = ~isnan(amp) & ~isnan(ampM) & ~isnan(choice) & ~isnan(allCoh);
            y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
            % x = [amp(utrls)' allCoh(utrls)' ampM(utrls)' subjMat(:, utrls)'];
            x = [amp(utrls)' allCoh(utrls)' ampM(utrls)'];
            
            [b, dev, stat] = glmfit(x,y,...
                'binomial', 'link', 'logit');
            
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess{dim, dimM}(tt, :) = stat.beta;
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_p_sess{dim, dimM}(tt, :) = stat.p;
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_T_sess{dim, dimM}(tt, :) = stat.t;
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_SE_sess{dim, dimM}(tt, :) = stat.se;
            dim_predChoPool.Cho_xMedRes_norm_coh_d1_VarEx_sess{dim, dimM}(tt, :) = 1-var(stat.resid)./var(y);
            
            % Choice - zscored SL: partial, cohMat
            utrls = ~isnan(amp) & ~isnan(ampM) & ~isnan(choice) & ~isnan(allCoh);
            y = choice(utrls)'; %y = Cho(~isnan(resAmp))';
            % x = [amp(utrls)' cohMat(:, utrls)' ampM(utrls)' subjMat(:, utrls)'];
            x = [amp(utrls)' cohMat(:, utrls)' ampM(utrls)'];
            
            [b, dev, stat] = glmfit(x,y,...
                'binomial', 'link', 'logit');
            
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_B_sess{dim, dimM}(tt, :) = stat.beta;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_p_sess{dim, dimM}(tt, :) = stat.p;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_T_sess{dim, dimM}(tt, :) = stat.t;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_SE_sess{dim, dimM}(tt, :) = stat.se;
            dim_predChoPool.Cho_xMedRes_norm_cohMat_d1_VarEx_sess{dim, dimM}(tt, :) = 1-var(stat.resid)./var(y);
        end
    end
end

% save('E:\Matalb analyses\dim_predChoPool', 'dim_predChoPool', 'par_sess')


figure; hold on; 
plot(whichTs, dim_predChoPool.ChoRes_normB_coh_d1_sess{dim}(:, 2), 'LineWidth', 2)
plot(whichTs, dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess{dim, 1}(:, 2))
plot(whichTs, dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess{dim, 2}(:, 2))
plot(whichTs, dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess{dim, 31}(:, 2))
plot(whichTs, dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess{dim, 32}(:, 2))
plot(whichTs, dim_predChoPool.Cho_xMedRes_norm_coh_d1_B_sess{dim, 20}(:, 2))
legend('raw', 'tin-c', 'tin-i', 'min-c', 'min-i', 'pc1')

% RT ----------------------------------------------------------------
medDims = useDims;
for cho = 1 : 2
    for dim = useDims
        for dimM = medDims
            disp([num2str(dim) ' vs ' num2str(dimM)])
            
            for tt = 1 : length(whichTs)
                amp = [];
                ampM = [];
                allCoh = [];
                rt = [];
                % session = [];
                for sess = 1 : 8
                    trls = find(resDim_norm_choice{sess} == 2-cho);
                    addAmp = resDimT_norm{dim, sess}(trls, tt)';
                    amp = [amp addAmp];
                    addAmp = resDimT_norm_compAmp{dimM, sess}(1, trls);
                    ampM = [ampM addAmp];
                    c = resDim_norm_sigcoh{sess}(1, trls);
                    allCoh = [allCoh c];
                    c = resDim_norm_rt{sess}(1, trls);
                    rt = [rt c];
                    % session = [session sess * ones(1, length(c))];
                end
                cohs = unique(allCoh);
                cohMat = [];
                for c = 1 : length(cohs)
                    cohMat(c, :) = (allCoh == cohs(c));
                end
                %             subjMat = [];
                %             for c = unique(session)
                %                 subjMat(c, :) = (session == c);
                %             end
                
                if dim == dimM
                    % RT - SL normalized (not zscored)
                    utrls = ~isnan(amp) & ~isnan(rt);
                    y = rt(utrls)';
                    x = amp(utrls)';
                    
                    [R,P] = corrcoef(x, y, 'rows', 'pairwise');
                    
                    dim_predXPool.pre_normRTr{dim, cho}(tt) = R(1, 2);
                    dim_predXPool.pre_normRTp{dim, cho}(tt) = P(1, 2);
                    dim_predXPool.pre_normRT_N(cho) = length(x);
                    
                    % RT - SL normalized with coherence indicator (not zscored)
                    utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt);
                    y = rt(utrls)';
                    x = amp(utrls)';
                    z = allCoh(utrls)';
                    
                    [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                    
                    dim_predXPool.pre_normCohRTr{dim, cho}(tt) = rho;
                    dim_predXPool.pre_normCohRTp{dim, cho}(tt) = p;
                    
                    % RT - SL normalized with coherence indicator (not zscored)
                    utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt);
                    y = rt(utrls)';
                    x = amp(utrls)';
                    z = cohMat(:, utrls)';
                    
                    [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                    
                    dim_predXPool.pre_normCohMatRTr{dim, cho}(tt) = rho;
                    dim_predXPool.pre_normCohMatRTp{dim, cho}(tt) = p;
                    
                    
                end
                % RT - SL normalized (not zscored) - mediated
                utrls = ~isnan(amp) & ~isnan(rt) & ~isnan(ampM);
                y = rt(utrls)';
                x = amp(utrls)';
                z = ampM(utrls)';
                
                
                [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                
                dim_predXPool.preX_normRTr{dim, dimM, cho}(tt) = rho;
                dim_predXPool.preX_normRTp{dim, dimM, cho}(tt) = p;
                
                % RT - SL normalized with coherence indicator (not zscored) - mediated
                utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt) & ~isnan(ampM);
                y = rt(utrls)';
                x = amp(utrls)';
                z = [allCoh(utrls)' ampM(utrls)'];
                
                [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                
                dim_predXPool.preX_normCohRTr{dim, dimM, cho}(tt) = rho;
                dim_predXPool.preX_normCohRTp{dim, dimM, cho}(tt) = p;
                
                % RT - SL normalized with coherence indicator (not zscored) - mediated
                utrls = ~isnan(amp) & ~isnan(allCoh) & ~isnan(rt) & ~isnan(ampM);
                y = rt(utrls)';
                x = amp(utrls)';
                z = [cohMat(:, utrls)' ampM(utrls)'];
                
                [rho, p] = partialcorr(x, y, z, 'rows', 'pairwise');
                
                dim_predXPool.preX_normCohMatRTr{dim, dimM, cho}(tt) = rho;
                dim_predXPool.preX_normCohMatRTp{dim, dimM, cho}(tt) = p;
                
            end
        end
    end
end

% save('E:\Matalb analyses\dim_pred_pool', 'dim_predChoPool', 'dim_predXPool', 'par_sess', '-v7.3')

figure; hold on; 
plot(whichTs, dim_predXPool.pre_normRTr{dim, cho},'r', 'LineWidth', 2)
plot(whichTs, dim_predXPool.pre_normCohRTr{dim, cho},'g', 'LineWidth', 2)
plot(whichTs, dim_predXPool.pre_normCohMatRTr{dim, cho},'b', 'LineWidth', 2)
plot(whichTs, dim_predXPool.preX_normRTr{dim, dimM, cho},'r', 'LineWidth', 1)
plot(whichTs, dim_predXPool.preX_normCohRTr{dim, dimM, cho},'g', 'LineWidth', 1)
plot(whichTs, dim_predXPool.preX_normCohMatRTr{dim, dimM, cho},'b', 'LineWidth', 1)
legend('raw', 'tin-c', 'tin-i', 'min-c', 'min-i', 'pc1')

% 
dim = 1; cho = 1;
figure; hold on; 
plot(whichTs, dim_predXPool.pre_normCohMatRTr{dim, cho},'k', 'LineWidth', 2)
dimM=1; plot(whichTs, dim_predXPool.preX_normCohMatRTr{dim, dimM, cho}, 'LineWidth', 1)
dimM=2; plot(whichTs, dim_predXPool.preX_normCohMatRTr{dim, dimM, cho}, 'LineWidth', 1)
dimM=31; plot(whichTs, dim_predXPool.preX_normCohMatRTr{dim, dimM, cho}, 'LineWidth', 1)
dimM=32; plot(whichTs, dim_predXPool.preX_normCohMatRTr{dim, dimM, cho}, 'LineWidth', 1)
dimM=20; plot(whichTs, dim_predXPool.preX_normCohMatRTr{dim, dimM, cho}, 'LineWidth', 1)
legend('raw', 'tin-c', 'tin-i', 'min-c', 'min-i', 'pc1')


%% Compute separately per session again
% Make sure that scaling works out across sessions --> z-score each time
% point

%% Compute correlation with RT


% -----------------------------------
% Setting used in analyses for paper
par = par_sess{1};
whichTs = [-0.2 : 0.005 : 0.6]; 
compBin = find(par.tsl >=0.525 & par.tsl <=0.575);
t_CB = [.525 .575]; t_CBLate = [.525 .575]; 
% Do I need the three below?
t_CBr = round([-.1 -0.05], 3); cbR = [find(round(par.trl, 3) == t_CBr(1)) find(round(par.trl, 3) ==t_CBr(2))];
t_CBresp = par.tsl(cbR);
t_BL =  [-.1 0];

rtAdd = 0.08; % used to be: rtAdd = 0.1;
bWidth = 0.08; minRT = [0.55+bWidth/2+rtAdd]; minRT_comp = minRT;
minRT_min = 0.475; 
maxRT = 2; %1.6;
useS = 1 : 10;

% -----------------------------------

whichTsPart = whichTs;
he = 1; 
useRes = 1; useZres = 1; useRegress = 1;

% sess = 1; useSC = 1; 
calcSess = [1 : 8];

dim_pred = [];
dim_pred.whichTs = whichTs;
dim_pred.compBin = compBin;
dim_pred.whichTsPart = whichTsPart;
dim_pred.bWidth = bWidth;
dim_pred.t_CBLate = t_CBLate;
dim_pred.rtAdd = rtAdd;
dim_pred.minRT_comp = minRT_comp;
dim_pred.minRT_min = minRT_min;
dim_pred.maxRT = maxRT;
dim_pred.t_CBr = t_CBr;
dim_pred.t_CBresp = t_CBresp;
dim_pred.t_BL = t_BL;


dim_pred.RTr = nan(length(calcSess), length(allDim_sl), 2, 4, 2, length(whichTs)); % Do I want to make this big just for the two dendro dimensions? 
dim_pred.RTp = dim_pred.RTr;
dim_pred.RTrPart = dim_pred.RTr;
dim_pred.RTpPart = dim_pred.RTr;
dim_pred.preRTrRes = nan(length(calcSess), length(allDim_sl), 4, 2, length(whichTs));
dim_pred.preRTpRes = dim_pred.preRTrRes;
dim_pred.preRTregressB = nan(length(calcSess), length(allDim_sl), 4, 2, length(whichTs), 3);
dim_pred.preRTregressP = dim_pred.preRTregressB;
dim_pred.preRTrPartCoh = dim_pred.preRTrRes;
dim_pred.preRTpPartCoh = dim_pred.preRTrRes;
dim_pred.preRTrPartCohBL = dim_pred.preRTrRes;
dim_pred.preRTpPartCohBL = dim_pred.preRTrRes;
dim_pred.preRTrResPart = dim_pred.preRTrRes;
dim_pred.preRTrResZPart = dim_pred.preRTrRes;
dim_pred.preRTrPartT2Coh = dim_pred.preRTrRes;
dim_pred.preRTpPartT2Coh = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohBL = dim_pred.preRTrRes;

dim_pred.RTrPartTin = dim_pred.RTr;
dim_pred.RTpPartTin = dim_pred.RTr;
dim_pred.preRTrResPartTin = dim_pred.preRTrRes;
dim_pred.preRTrResZPartTin = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohTin = dim_pred.preRTrRes;
dim_pred.preRTpPartT2CohTin = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohBLTin = dim_pred.preRTrRes;
dim_pred.preRTpPartT2CohBLTin = dim_pred.preRTrRes;

dim_pred.RTrPartTinT = dim_pred.RTr;
dim_pred.RTpPartTinT = dim_pred.RTr;
dim_pred.preRTrResPartTinT = dim_pred.preRTrRes;
dim_pred.preRTrResZPartTinT = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohTinT = dim_pred.preRTrRes;
dim_pred.preRTpPartT2CohTinT = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohBLTinT = dim_pred.preRTrRes;
dim_pred.preRTpPartT2CohBLTinT = dim_pred.preRTrRes;

dim_pred.RTrPartresp = dim_pred.RTr;
dim_pred.RTpPartresp = dim_pred.RTr;
dim_pred.preRTrResPartresp = dim_pred.preRTrRes;
dim_pred.preRTrResZPartresp = dim_pred.preRTrRes;
dim_pred.preRTrPartT2Cohresp = dim_pred.preRTrRes;
dim_pred.preRTpPartT2Cohresp = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohBLresp = dim_pred.preRTrRes;

dim_pred.RTrPartToutT = dim_pred.RTr;
dim_pred.RTpPartToutT = dim_pred.RTr;
dim_pred.preRTrResPartToutT = dim_pred.preRTrRes;
dim_pred.preRTrResZPartToutT = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohToutT = dim_pred.preRTrRes;
dim_pred.preRTpPartT2CohToutT = dim_pred.preRTrRes;
dim_pred.preRTrPartT2CohBLToutT = dim_pred.preRTrRes;
dim_pred.preRTpPartT2CohBLToutT = dim_pred.preRTrRes;


dim_pred.RTr_noMin = dim_pred.RTr;
dim_pred.RTp_noMin = dim_pred.RTr;
dim_pred.preRTrRes_noMin = nan(10, length(allDim_sl), 4, 2, length(whichTs));
dim_pred.preRTpRes_noMin = nan(10, length(allDim_sl), 4, 2, length(whichTs));
dim_pred.preRTrPartCohBL_noMin = dim_pred.preRTrRes;
dim_pred.preRTpPartCohBL_noMin = dim_pred.preRTrRes;

% Redistribute parameters
whichTs = dim_pred.whichTs;
compBin = dim_pred.compBin;
whichTsPart = dim_pred.whichTsPart;
bWidth = dim_pred.bWidth;
t_CBLate = dim_pred.t_CBLate;
rtAdd = dim_pred.rtAdd;
minRT_comp = dim_pred.minRT_comp;
minRT_min = dim_pred.minRT_min;
maxRT = dim_pred.maxRT; % 2
t_CBr = dim_pred.t_CBr;
t_CBresp = dim_pred.t_CBresp;
t_BL = dim_pred.t_BL;

he = 1; useRes = 1; useZres = 0; useRegress = 1;

calcDims = [1 24 25 50];
calcDims = [1 2 3 4 40 41];
calcDims = [20];
calcDims = [24];
calcDims = [57 58];
calcDims = [26];
calcDims = [1 2 3 4 20 24 26 40 41 57 58];
calcDims = [26 40 41 57 58];
calcDims = [1 2 3 4 20 24];
calcDims = [59 60 61];
calcDims = [24];
calcDims = [1 2 3 20 24 26 41 59 61];
calcDims = 59;
calcDims = [31 32 33];
calcDims = [62 63];
calcDims = [21];
calcDims = [64 65]; % forgot to save these dimensions!
calcDims = [66]; % forgot to save these dimensions!
calcDims = [67];
calcDims = [68 69];

calcDims = [1 2 7 20 31 32 61 63];


for sess = calcSess
    for dim = calcDims %24 : 25 %[1 : length(allDim_sl)]
        disp(num2str(dim))
        if size(allDim_sl{dim}, 1) >= sess | size(allDim_sl{dim}, 2)
            par = par_sess{sess};
            dataMat = par_sess{sess}.dataMat;
            chos = unique(dataMat(:, par.idx.cho_trg))';

            if dim >= 16 & dim <=23 % PCA
                useSC = 1 : 3; %pca_intDim{dim}(sess);
            elseif dim == 8 | dim == 10 
                useSC = 1 : size(allDim_sl{dim}, 2);
            elseif dim == 24 | dim == 25
                useSC = 1 : 4; % size(allDim_sl{dim}, 1);
                elseif dim == 61
                useSC = 1 : 3; % size(allDim_sl{dim}, 1);
                elseif dim == 59 | dim == 60
                useSC = 1 : 5; % size(allDim_sl{dim}, 1);
            else
                useSC = 1;
            end
            dim_pred.useSC{dim, sess} = useSC;
            
            rawFRs = []; rawFRsCompS = []; rawFRsCompS_Tout = [];
            for sc = 1 : length(useSC)
                if length(useSC) > 1
                    if dim == 8 | dim == 10
                        rawFRs{sc} = allDim_sl{dim}{sess, useSC(sc)};
                    elseif dim >= 16 & dim <=23 % PCA
                        rawFRs{sc} = allDim_sl{dim}{useSC(sc), sess};
                    elseif dim == 59 | dim == 60 % dPCA
                        rawFRs{sc} = allDim_sl{dim}{useSC(sc), sess};
                    else
                        rawFRs{sc} = allDim_sl{dim}{useSC(sc), sess};
                    end
                else
                    rawFRs{sc} = allDim_sl{dim}{useSC(sc), sess};
                end
                rawFRsCompS{sc} = allDim_sl{1}{sess}; t_CB = t_CBLate;
                rawFRsCompS_Tout{sc} = allDim_sl{2}{sess}; t_CB = t_CBLate;
                % rawFRsCompS_DD{sc} = allDim_sl{3}{sess}; t_CB = t_CBLate;
                % rawFRsCompS_ED{sc} = allDim_sl{5}{sess}; t_CB = t_CBLate;
                
            end
            if sum(sum(isnan(rawFRs{sc}))) ~= length(rawFRs{sc}(:))
                for justZero = 0 : 1
                    ccc = 1 : length(useSC);
                    % Shorter minimum RT
                    minRT = minRT_min;
                    calcPredRT_earlyAmp_clean
                    
                    dim_pred.RTr_noMin(sess, dim, justZero+1, ccc, :, :) = preRTr;
                    dim_pred.RTp_noMin(sess, dim, justZero+1, ccc, :, :) = preRTp;
                    
                    if justZero == 0
                        dim_pred.preRTrPartCohBL_noMin(sess, dim, ccc, :, :) = preRTrPartCohBL;
                        dim_pred.preRTpPartCohBL_noMin(sess, dim, ccc, :, :) = preRTpPartCohBL;
                        
                        if useRes == 1
                            dim_pred.preRTrRes_noMin(sess, dim, ccc, :, :) = preRTrRes;
                            dim_pred.preRTpRes_noMin(sess, dim, ccc, :, :) = preRTpRes;
                            
                            dim_pred.preRTrResZ_noMin(sess, dim, ccc, :, :) = preRTrResZ;
                            dim_pred.preRTpResZ_noMin(sess, dim, ccc, :, :) = preRTpResZ;
                        end
                        
                    end
                    
                    % Minimim RT = 0.65s
                    minRT = minRT_comp;
                    calcPredRT_earlyAmp_clean;
                    
                    dim_pred.RTr(sess, dim, justZero+1, ccc, :, :) = preRTr;
                    dim_pred.RTp(sess, dim, justZero+1, ccc, :, :) = preRTp;
                    
                    if justZero == 0
                        if useRes == 1
                            dim_pred.preRTrRes(sess, dim, ccc, :, :) = preRTrRes;
                            dim_pred.preRTpRes(sess, dim, ccc, :, :) = preRTpRes;
                            
                            dim_pred.preRTregressB(sess, dim, ccc, :, :, :) = preRTregressB;
                            dim_pred.preRTregressP(sess, dim, ccc, :, :, :) = preRTregressP;
                            
                            dim_pred.preRTrResZ(sess, dim, ccc, :, :) = preRTrResZ;
                            dim_pred.preRTpResZ(sess, dim, ccc, :, :) = preRTpResZ;
                        end
                        dim_pred.preRTrPartCoh(sess, dim, ccc, :, :) = preRTrPartCoh;
                        dim_pred.preRTpPartCoh(sess, dim, ccc, :, :) = preRTpPartCoh;
                        
                        dim_pred.preRTrPartCohBL(sess, dim, ccc, :, :) = preRTrPartCohBL;
                        dim_pred.preRTpPartCohBL(sess, dim, ccc, :, :) = preRTpPartCohBL;
                        
                    end

                    disp([num2str(sess) ': partial PredRT_earlyAmp self late...'])
                    t_CB = t_CBLate;
                    calcPredRT_earlyAmp_Partial_clean
                    
                    dim_pred.RTrPart(sess, dim, justZero+1, ccc, :, :) = preRTrPart;
                    
                    if justZero == 0
                        dim_pred.preRTrResPart(sess, dim, ccc, :, :) = preRTrResPart; % This one!
                        dim_pred.preRTpResPart(sess, dim, ccc, :, :) = preRTpResPart;
                        dim_pred.preRTrResZPart(sess, dim, ccc, :, :) = preRTrResZPart;
                        dim_pred.preRTrPartT2Coh(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTpPartT2Coh(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTrPartT2CohBL(sess, dim, ccc, :, :) = preRTrPartT2CohBL;
                        dim_pred.preRTpPartT2CohBL(sess, dim, ccc, :, :) = preRTpPartT2CohBL;
                    end
                    
                    disp([num2str(sess) ': Calculating partial PredRT_earlyAmp Tin late...'])
                    latePartBin = 1; rawFRsComp = rawFRsCompS;
                    calcPredRT_earlyAmp_Partial_cross
                    
                    dim_pred.RTrPartTin(sess, dim, justZero+1, ccc, :, :) = preRTrPart;
                    dim_pred.RTpPartTin(sess, dim, justZero+1, ccc, :, :) = preRTpPart;
                    
                    if justZero == 0
                        if useRes == 1
                            dim_pred.preRTrResPartTin(sess, dim, ccc, :, :) = preRTrResPart; % This one!
                            dim_pred.preRTpResPartTin(sess, dim, ccc, :, :) = preRTpResPart;
                            dim_pred.preRTrResZPartTin(sess, dim, ccc, :, :) = preRTrResZPart;
                            dim_pred.preRTpResZPartTin(sess, dim, ccc, :, :) = preRTpResZPart;
                        end
                        dim_pred.preRTrPartT2CohTin(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTpPartT2CohTin(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTrPartT2CohBLTin(sess, dim, ccc, :, :) = preRTrPartT2CohBL;
                        dim_pred.preRTpPartT2CohBLTin(sess, dim, ccc, :, :) = preRTpPartT2CohBL;
                    end
                    
                    
                    disp([num2str(sess) ': partial PredRT_earlyAmp Tin same...'])
                    latePartBin = 0;
                    calcPredRT_earlyAmp_Partial_cross
                    
                    dim_pred.RTrPartTinT(sess, dim, justZero+1, ccc, :, :) = preRTrPart;
                    
                    if justZero == 0
                        if useRes == 1
                            dim_pred.preRTrResPartTinT(sess, dim, ccc, :, :) = preRTrResPart; % This one?
                            dim_pred.preRTpResPartTinT(sess, dim, ccc, :, :) = preRTpResPart;
                            dim_pred.preRTrResZPartTinT(sess, dim, ccc, :, :) = preRTrResZPart;
                            dim_pred.preRTpResZPartTinT(sess, dim, ccc, :, :) = preRTpResZPart;
                        end
                        dim_pred.preRTrPartT2CohTinT(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTpPartT2CohTinT(sess, dim, ccc, :, :) = preRTpPartT2Coh;
                        dim_pred.preRTrPartT2CohBLTinT(sess, dim, ccc, :, :) = preRTrPartT2CohBL;
                        dim_pred.preRTpPartT2CohBLTinT(sess, dim, ccc, :, :) = preRTpPartT2CohBL;
                    end
                    
                    disp([num2str(sess) ': partial PredRT_earlyAmp Tout late...'])
                    latePartBin = 1; rawFRsComp = rawFRsCompS_Tout;
                    calcPredRT_earlyAmp_Partial_cross
                    
                    dim_pred.RTrPartTout(sess, dim, justZero+1, ccc, :, :) = preRTrPart;
                    
                    if justZero == 0
                        if useRes == 1
                            dim_pred.preRTrResPartTout(sess, dim, ccc, :, :) = preRTrResPart;
                            dim_pred.preRTrResZPartTout(sess, dim, ccc, :, :) = preRTrResZPart;
                        end
                        dim_pred.preRTrPartT2CohTout(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTpPartT2CohTout(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTrPartT2CohBLTout(sess, dim, ccc, :, :) = preRTrPartT2CohBL;
                        dim_pred.preRTpPartT2CohBLTout(sess, dim, ccc, :, :) = preRTpPartT2CohBL;
                    end
                    
                    disp([num2str(sess) ': partial PredRT_earlyAmp Tout same...'])
                    latePartBin = 0; rawFRsComp = rawFRsCompS_Tout;
                    calcPredRT_earlyAmp_Partial_cross
                    
                    dim_pred.RTrPartToutT(sess, dim, justZero+1, ccc, :, :) = preRTrPart;
                    
                    if justZero == 0
                        if useRes == 1
                            dim_pred.preRTrResPartToutT(sess, dim, ccc, :, :) = preRTrResPart;
                            dim_pred.preRTrResZPartToutT(sess, dim, ccc, :, :) = preRTrResZPart;
                        end
                        dim_pred.preRTrPartT2CohToutT(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTpPartT2CohToutT(sess, dim, ccc, :, :) = preRTrPartT2Coh;
                        dim_pred.preRTrPartT2CohBLToutT(sess, dim, ccc, :, :) = preRTrPartT2CohBL;
                        dim_pred.preRTpPartT2CohBLToutT(sess, dim, ccc, :, :) = preRTpPartT2CohBL;
                    end
                    
                    clear preRTr preRTp preRTrRes preRTpRes  preRTregressB preRTregressP preRTrPartCoh preRTpPartCoh preRTrPartCohBL preRTpPartCohBL
                    clear preRTrPart preRTrResPart preRTrResZPart preRTrPartT2Coh preRTrPartT2Coh preRTrPartT2CohBL
                    clear preRTpPart preRTpResPart preRTpResZPart preRTpPartT2Coh preRTpPartT2Coh preRTpPartT2CohBL preRTpResZ preRTpReszPart
                end
            end
        end
    end
end

% note = ['S1-8. Just dims [1 : 5 20 24 26 31 32 33 40 41 57 58 59 60 : 66]. ',...
%     'all dims updated 220423: partials did not include max RT, will have to recalculate all...',... 
%     ' redid 59 with updated dPCA 220426.',...
%     'Last saved 230103 up including 64:69.']; 
% save('E:\Matalb analyses\sl_dim_pred_TinDD_S1-10', 'dat_sess', 'par_sess',...
%     'dim_pred', 'allW','dimName', 'note', '-v7.3') % par_sess is updated!



%% Choice




%% Plot RT across session mean and SEM

% load('E:\Matalb analyses\RTvAmp_xDim_230102') 
% % 'dat_sess', 'par_sess','preRTrPartT2CohBL_cross', 'preRTpPartT2CohBL_cross', 'preRTrPartT2CohBL_crossT', 'preRTpPartT2CohBL_crossT','preRTrResPart_cross', 'preRTpResPart_cross'
% load('E:\Matalb analyses\sl_dim_pred_TinDD_S1-10')

whichTsPart = dim_pred.whichTsPart;
minRT_comp = dim_pred.minRT_comp;
minRT = minRT_comp;
bWidth = dim_pred.bWidth;
maxRT = dim_pred.maxRT;
whichTs = dim_pred.whichTs;
compBin = dim_pred.compBin;
rtAdd = dim_pred.rtAdd;
t_CBLate = dim_pred.t_CBLate;
t_CBr = dim_pred.t_CBr;
t_CBresp = dim_pred.t_CBresp;
t_BL = dim_pred.t_BL;
pltDim = dimName; pltDim{1} = 'Tin-c'; pltDim{2} = 'Tin-i'; pltDim{20} = 'PC1'; pltDim{63} = 'rampD'; pltDim{69} = 'Tin w_r_a_m_p';

plotTraw = find(whichTs >=0.2 & whichTs <=0.45);
plotTblock = plotTraw;

xLims = [0.2 0.5]; yLims = [-0.35 0.05];
compTP = find(round(whichTs, 3) == 0.55);
plotSess = [1 : 8]; % [1 3 : 8]; % 
barCols = {[233 149 53]/255, [76 0 153]/255, [0 153 76]/255, [233 53 149]/255,...
    [0 149 233]/255, [76 0 153]/255, [0 153 76]/255, [149 53 149]/255};

plotDims = [1 2 20 63 69]; medDims = plotDims;
plotDims = [1 2 31 32 20 63]; medDims = [1 2 31 32]; 
plotDims = [63 20 1]; medDims = [63 20 1]; 
yLims = [-0.4 0.1]; yLims = [-0.35 0.1];
rows = 3;
colms = ceil(length(plotDims)/rows);

lw = 2;
faceAlpha = 1;

figure('Position', [200 10 1500 1000]); hold on
for d1 = 1 : length(plotDims)
    subplot(rows, colms, d1); hold on
    dim = plotDims(d1); title(pltDim{dim})
    cho = 1; if dim == 57 | dim == 58; cho = 2; end
    
    rtR = nan(8, length(whichTs));
    for sess = plotSess
        sc = 1;
        % if dim == 20; sc = 1; end
        rtR(sess, :) = squeeze(dim_pred.preRTrRes(sess, dim, sc, cho, :));
    end
    pltMeanRaw = nanmean(rtR, 1);
    pltSEMRaw = nanstd(rtR, [], 1) / sqrt(size(rtR, 1));
    rtMraw = nanmean(rtR(:, compTP), 1);

    %     rtR = nan(8, length(whichTs));
    %     for sess = plotSess
    %         rtR(sess, :) = squeeze(dim_pred.preRTrResPart(sess, dim, sc, cho, :));
    %     end
    %
    %     pltMeanBlk = nanmean(rtR, 1);
    %     pltSEMBlk = nanstd(rtR, [], 1) / sqrt(size(rtR, 1));
    
    plotCol = 'k';
    clear s
    s = shadedErrorBar(whichTs(plotTraw), pltMeanRaw(plotTraw), pltSEMRaw(plotTraw), 'transparent', faceAlpha); %'facealpha',faceAlpha);
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lw;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    %plot(whichTs(compTP), rtMraw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol, 'MarkerSize', 12)
    %rtSEM = nanstd(rtR(:, compTP), [], 1) / sqrt(size(rtR(:, compTP), 1));
    %line(whichTs(compTP) * [1 1], rtMraw + [-rtSEM rtSEM], 'Color', plotCol)
    
    %         plotCol = [255 69 0]/255; %[250 128 114]/255; % [128 0 128]/255; %
    %         clear s
    %         s = shadedErrorBar(whichTs(plotTblock), pltMeanBlk(plotTblock), pltSEMBlk(plotTblock));
    %         set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    %         s.mainLine.LineWidth = lw;
    %         s.patch.FaceColor = plotCol;
    %         s.mainLine.Color = plotCol;
    
    for d2 = 1 : length(medDims)
        dimM = medDims(d2);
        clear s
        if dim ~= dimM
            rtR = nan(8, length(whichTs));
            for sess = plotSess
                rtR(sess, :) = preRTrResPart_cross{dim, dimM, sess}(1, cho, :);
            end
        else
            rtR = nan(8, length(whichTs));
            for sess = plotSess
                rtR(sess, :) = squeeze(dim_pred.preRTrResPart(sess, dim, sc, cho, :));
            end
        end
        plotCol = barCols{d2};
        pltMean = nanmean(rtR, 1);
        pltSEM = nanstd(rtR, [], 1) / sqrt(size(rtR, 1));
        s = shadedErrorBar(whichTs(plotTblock), pltMean(plotTblock), pltSEM(plotTblock), 'transparent', faceAlpha); % 'facealpha',faceAlpha);
        % plot(whichTs(plotTblock), pltMean(plotTblock), 'Color', plotCol, 'LineWidth', 1.5)
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = lw;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        %             pltTinComp = squeeze(dim_pred.preRTrRes(:, 1, sc, cho, compTP));
        %             rtM = nanmean(rtR(:, compTP), 1);
        %             plot(whichTs(compTP), rtM, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol, 'MarkerSize', 12)
        %             rtSEM = nanstd(rtR(:, compTP), [], 1) / sqrt(size(rtR(:, compTP), 1));
        %             line(whichTs(compTP) * [1 1], rtM + [-rtSEM rtSEM], 'Color', plotCol)
        
        
        
    end
    ylabel('Correlation (r)')
    xlabel('Time after motion onset [s]')
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLims); if dim == 2; ylim(-flip(yLims)); end
    set(gca,'TickDir','out', 'FontSize', 10);
    
end
% legend('Raw', 'Tin-c', 'Tin-i', 'PC1', 'rampD', 'ramp Tin only')
% legend('Raw', 'Tin-c', 'Tin-i', 'Min-c', 'Min-i')
legend('Raw', 'Ramp', 'PC1', 'TinC')







