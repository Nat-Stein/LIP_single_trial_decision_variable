% calc_PCA_Min

%% ----------------------------------------
% PCA analyses on LIP single-trial data
%% ----------------------------------------


% Remove sl_pca and rl_pca?


%% Load data
load(fullfile(saveLoc, 'par_sess'))

%% Generate across trial averages and repeat PCA

corr = 1; addRL = 0;
preCut = 100; tic
par = par_sess{1};
par.tsl = round(par.tsl,3);
useT = [0.2 0.6]; bDist = 50; bWidth = 50; bWidthI = 1;
t_incl = find(par.tsl >= useT(1) & par.tsl <= useT(2));
tpts = find(par.tsl==useT(1)) : bWidth : find(par.tsl==useT(2));
thisT = round(par.tsl, 3);

for sess = 1 : 8
    
    % ----------------------------------------------
    % Load session specific data
    disp(['Session ' num2str(sess) ': Loading data...'])
    % Use data normalized based on mean and std in one time window
    % 200-600ms after motion onset
    suff = '_mean_stp50';
    clear spkZnan
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZnan')
    disp(['Session ' num2str(sess) ': Pre-processing...'])
    par = par_sess{sess};
    goodN = [par.unit_MinC par.unit_MinI];
    dataMat_SL = spkZnan.SL(goodN, :, :);
    dataMat_RL = spkZnan.RL(goodN, :, :);
    
    % ----------------------------------------------
    % Assign behavioral parameters for specific session
    
    beh = par.beh;
    par.tsl = round(par.tsl, 3);
    
    % ----------------------------------------------
    % Compute across trial averages within signed coherence
    allConds = nan(12, 2, size(dataMat_SL, 1), size(dataMat_SL, 3));
    allCconcat = []; allCconcatRL = []; conds = []; cc = 0;
    
    sigCohs = unique(beh.sig_coh);
    sigCohs = [sigCohs 0]; % Adding a second "0" to be able to split 0%-coherence trials with left and right choices
    zeroS = 0;
    for coh = 1 : length(sigCohs)
        if sigCohs(coh) == 0
            % Split 0% coherenec trials by choice 
            trls = find(beh.sig_coh == sigCohs(coh) & beh.cho_trg == zeroS);
            zeroS = zeroS + 1;
        else
            % For coherences other than 0%, include only correct trials
            trls = find(beh.sig_coh == sigCohs(coh) & beh.correct == corr);
            % trls = find(beh.sig_coh == sigCohs(coh));
        end
        nTRL(coh, corr+1) = length(trls);
        if length(trls) > 5
            allConds(coh, corr+1, :, :) = squeeze(nanmean(dataMat_SL(:, trls, :), 2));
            cc = cc + 1; allCconcat(cc, :, :) = squeeze(nanmean(dataMat_SL(:, trls, :), 2));
            conds(cc, :) = [coh corr];
            medianRTs(cc) = round(nanmedian(beh.rt(trls)), 3);
        end
    end
    
    thisSLmat_M = permute(allCconcat, [2, 1, 3]);
    
    % ----------------------------------------------
    % Save input data
    inMean = [];
    inMean.allConds = allConds;
    inMean.allCconcat = allCconcat;
    inMean.conds = conds;
    inMean.medianRTs = medianRTs;
    
    % ----------------------------------------------
    % Concatenate across-trial averages for all coherences
    dataXdot = [];
    nTrialD = [];
    nCoh = [];
    for c = 1 : size(thisSLmat_M, 2)
        thisRT = find(round(thisT, 3) == round(medianRTs(c), 3));
        for tpt =  1 : length(tpts)
            % Define the end of the time bin to exclude any data in the
            % last 100ms before the saccade and any time beyond 600ms
            maxT = min(round(thisRT - preCut - (bWidth*bWidthI)/2), round(tpts(tpt) + bWidth/2));
            
            tVec = [round(tpts(tpt) - bWidth/2) : maxT];
            if length(tVec) >= bWidth/5
                % Delay period
                thisTRL = squeeze(nanmean(thisSLmat_M(:, c, tVec), 3));
                if sum(isnan(thisTRL)) == length(thisTRL)
                    break;
                end
                dataXdot = cat(2, dataXdot, thisTRL);
                
                nTrialD = [nTrialD c];
                nCoh = [nCoh conds(c, 1)];
                
            end
        end
    end
    
    % ----------------------------------------------
    % Compute PCA on RDM mean per condition
    disp(['Session ' num2str(sess) ': Computing PCA...'])
    [coeffD,scoreD,latentD,tsquaredD,explainedD,muD] = pca(dataXdot');
    
    % Compute single-trial DV for first x PCs
    disp(['Session ' num2str(sess) ': Computing projections...'])
    cLen = 100;
    for pc = 1 : min(10, size(coeffD,2))
        w = coeffD(:, pc);
        % invert PC1 for session 4 to match polarity of other directions (activity on contra trials is greater than on ipsi trials)
        % if sess == 4 && pc == 1; w = -w; end     
        wAll = zeros(1, size(spkZnan.SL, 1)); wAll(goodN) = w;
        Csl = nan(size(dataMat_SL, 2), size(dataMat_SL, 3));
        Crl = Csl;
        for cc = 1 : ceil(size(dataMat_SL, 3)/cLen)
            tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(dataMat_SL, 3))];
            Csl(:, tps) = squeeze(mmx('mult', w', dataMat_SL(:, :, tps)));
            Crl(:, tps) = squeeze(mmx('mult', w', dataMat_RL(:, :, tps)));
        end
        sl_pca{pc, sess} = Csl;
        rl_pca{pc, sess} = Crl;
        w_pca{pc, sess} = wAll;
    end
    outPCA.coeffD = coeffD;
    outPCA.scoreD = scoreD;
    outPCA.latentD = latentD;
    outPCA.tsquaredD = tsquaredD;
    outPCA.explainedD = explainedD;
    outPCA.muD = muD;
    outPCA.suff = suff;
    outPCA.useT = useT;
    outPCA.bDist = bDist;
    outPCA.bWidth = bWidth;
    outPCA.bWidthI = bWidthI;
    outPCA.preCut = preCut;
    outPCA.goodN = goodN;
    
    % ----------------------------------------------
    % Save PCA output for individual sessions including
    % sl_pca - Projection of data onto first 10 PCs stimulus-aligned
    % rl_pca - Projection of data onto first 10 PCs response-aligned
    % w_pca  - weight vectors for first 10 PCs
    % outPCA - output of PCA
    % inMean - input to PCA
    disp(['Session ' num2str(sess) ': Saving...'])
    save(fullfile(saveLoc, ['pca_Min_' suff '_S' num2str(sess)]), 'sl_pca', 'rl_pca', 'w_pca', 'outPCA', 'inMean', '-v7.3')

    w = w_pca{1, sess};
   save(fullfile(saveLoc, ['PC1_Min_w_S' num2str(sess)]), 'w')
end

a = [];
b = [];
for sess = 1 : 8
    a = cat(1, a, sl_pca{1, sess});
    b = cat(1, b, rl_pca{1, sess});
end
sigSL.PC1min = a;
sigRL.PC1min = b;

a = [];
b = [];
for sess = 1 : 8
    a = cat(1, a, sl_pca{2, sess});
    b = cat(1, b, rl_pca{2, sess});
end
sigSL.PC2min = a;
sigRL.PC2min = b;

figure; hold on
for sess = 1 : 8
    subplot(2, 4, sess); hold on
    tr = find(sigSL.sig_coh > 0.5 & sigSL.session == sess);
    plt = conv(nanmean(sigSL.PC1min(tr, :), 1), ones(1, 50), 'same');
    plt(par.tsl >= median(sigSL.rt(tr)) - 0.1) = nan;
    plot(par.tsl, plt)
    plt = conv(nanmean(sigSL.PC2min(tr, :), 1), ones(1, 50), 'same');
    plt(par.tsl >= median(sigSL.rt(tr)) - 0.1) = nan;
    plot(par.tsl, plt, '--')
    
    tr = find(sigSL.sig_coh < -0.5 & sigSL.session == sess);
    plt = conv(nanmean(sigSL.PC1min(tr, :), 1), ones(1, 50), 'same');
    plt(par.tsl >= median(sigSL.rt(tr)) - 0.1) = nan;
    plot(par.tsl, plt)
    plt = conv(nanmean(sigSL.PC2min(tr, :), 1), ones(1, 50), 'same');
    plt(par.tsl >= median(sigSL.rt(tr)) - 0.1) = nan;
    plot(par.tsl, plt, '--')
    xlim([-0.1 0.5])
end
legend('PC1 +51.2%', 'PC2 +51.2%', 'PC1 -51.2%', 'PC2 -51.2%')


%% Linear regression model to predict 


for sess = 1 : 8
    
    % ----------------------------------------------
    % Load session specific data
    disp(['Session ' num2str(sess) ': Loading data...'])
    % Use data normalized based on mean and std in one time window
    % 200-600ms after motion onset
    suff = '_mean_stp50';
    clear spkZnan
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZnan')
    disp(['Session ' num2str(sess) ': Pre-processing...'])
    par = par_sess{sess};
    goodN = [par.unit_MinC par.unit_MinI];
    dataMat_SL = spkZnan.SL(goodN, :, :);
    dataMat_RL = spkZnan.RL(goodN, :, :);
    
    % ----------------------------------------------
    % Assign behavioral parameters for specific session
    
    beh = par.beh;
    sig_coh = par.beh.sig_coh;
    rt = par.beh.rt;
    par.tsl = round(par.tsl, 3);
    numTr = size(dataMat_SL, 2);
    t_start = 0.1;
    
    % ----------------------------------------------
    % Set up matrix for regression by concatinating trials
    
    allCoh = [];
    allMin = [];
    for trl = 1 : numTr
        
        useT = find(par.tsl >= t_start & par.tsl <= rt(trl) - preCut/1000);
        
        allCoh = cat(2, allCoh, sig_coh(trl) * ones(1, length(useT)));
        
        newMin = squeeze(dataMat_SL(:, trl, useT));
        allMin = cat(2, allMin, newMin);
        
    end
    
    % Normalize input to regression
    
    allMin_m = nanmean(allMin, 2);
    allMin_s = nanstd(allMin')';
    
    allMin_norm = (allMin - repmat(allMin_m, 1, size(allMin, 2))) ./repmat(allMin_s, 1, size(allMin, 2));
    
    mdl = fitglm(allMin_norm',allCoh');
    
    w_out = table2array(mdl.Coefficients(2:end, 1));
    
    dataMat_SL_norm = (dataMat_SL - repmat(allMin_m, 1, numTr, size(dataMat_SL, 3))) ./ repmat(allMin_s, 1, numTr, size(dataMat_SL, 3));
    dataMat_RL_norm = (dataMat_RL - repmat(allMin_m, 1, numTr, size(dataMat_SL, 3))) ./ repmat(allMin_s, 1, numTr, size(dataMat_SL, 3));
    Csl = nan(size(dataMat_SL_norm, 2), size(dataMat_SL_norm, 3));
    Crl = Csl;
    cLen = 100;
    for cc = 1 : ceil(size(dataMat_SL_norm, 3)/cLen)
        tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(dataMat_SL_norm, 3))];
        Csl(:, tps) = squeeze(mmx('mult', w_out', dataMat_SL_norm(:, :, tps)));
        Crl(:, tps) = squeeze(mmx('mult', w_out', dataMat_RL_norm(:, :, tps)));
    end
    
    sl_minReg{sess} = Csl;
    rl_minReg{sess} = Crl;
    
    wAll = zeros(1, size(spkZnan.SL, 1)); wAll(goodN) = w_out;
    w_minReg{sess} = wAll;
    
    int_minReg(sess) = table2array(mdl.Coefficients(1, 1));
    
end




a = [];
b = [];
for sess = 1 : 8
    a = cat(1, a, sl_minReg{sess});
    b = cat(1, b, rl_minReg{sess});
end
sigSL.minReg = -a;
sigRL.minReg = -b;

figure; hold on
for sess = 1 : 8
    subplot(2, 4, sess); hold on
    tr = find(sigSL.sig_coh > 0.5 & sigSL.session == sess);
    plt = conv(nanmean(sigSL.minReg(tr, :), 1), ones(1, 50), 'same');
    plt(par.tsl >= median(sigSL.rt(tr)) - 0.1) = nan;
    plot(par.tsl, plt)
    
    tr = find(sigSL.sig_coh < -0.5 & sigSL.session == sess);
    plt = conv(nanmean(sigSL.minReg(tr, :), 1), ones(1, 50), 'same');
    plt(par.tsl >= median(sigSL.rt(tr)) - 0.1) = nan;
    plot(par.tsl, plt)
    xlim([-0.1 0.5])
end
% legend('PC1 +51.2%', 'PC2 +51.2%', 'PC1 -51.2%', 'PC2 -51.2%')

wMinC = [];
wMinI = [];
for sess = 1 : 8
    par = par_sess{sess};
    
    wMinC = [wMinC w_minReg{sess}(par.unit_MinC)];
    wMinI = [wMinI w_minReg{sess}(par.unit_MinI)];
    
end

sum(sign(wMinC)==-1) / length(wMinC)

sum(sign(wMinI)==1) / length(wMinI)



