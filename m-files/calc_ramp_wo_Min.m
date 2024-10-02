
for sess = 1 : 8
    
    disp([dNames{d} ', Session ' num2str(sess)])
    
    clear w
    ff = fullfile(saveLoc,...
        [dNames{d} '_w_S' num2str(sess)]);
    load(ff, 'w')
    if size(w, 1) > size(w, 2); w = w'; end
    % Take out Min neurons
    mins = [par_sess{sess}.unit_MinC par_sess{sess}.unit_MinI];
    
    w(mins) = 0;
    
    clear par dd beh behMap behView spk spkMap spkView media stdev
    clear spkZ spkMapZ spkViewZ spkZnan spkMapZnan spkViewZnan
    clear rl_rdm sl_rdm
    cLen = 200;
    
    %% Use z-scored data (ramp)  ---------------------------------
    disp(['Session ' num2str(sess) ': loading...'])
    load(fullfile(saveLoc, ['spk_S' num2str(sess)]))
    
    % Normalize data using saved mean and std before projecting
    load(fullfile(saveLoc, ['Ramp_w_S' num2str(sess)]), 'use_lasso_flag', 'media', 'stdev'); % Not loading: 'w','b0',
    
    goodN = find(stdev ~= 0);
    
    % Compute: xn = (x - media)./stdev;
    mMat = repmat(media(goodN)', 1, size(spk.SLnan, 2), size(spk.SLnan, 3));
    sMat = repmat(stdev(goodN)', 1, size(spk.SLnan, 2), size(spk.SLnan, 3));
    Cin = (squeeze(spk.SLnan(goodN, :, :)) - mMat) ./ sMat;
    
    disp(['Session ' num2str(sess) ': projecting data...'])
    % RDM: Stimulus-locked data
    C = nan(size(Cin, 2), size(Cin, 3));
    for cc = 1 : ceil(size(Cin, 3)/cLen)
        tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
        C(:, tps) = squeeze(mmx('mult', w(goodN), Cin(:, :, tps)));
    end
    sl_rdm = C; % figure; imagesc(C)
    % RDM: Response-locked data
    Cin = (squeeze(spk.RLnan(goodN, :, :)) - mMat) ./ sMat;
    C = nan(size(Cin, 2), size(Cin, 3));
    for cc = 1 : ceil(size(Cin, 3)/cLen)
        tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
        C(:, tps) = squeeze(mmx('mult', w(goodN), Cin(:, :, tps)));
    end
    rl_rdm = C;
    
    
    %% Save single session data matrices
    disp(['Session ' num2str(sess) ': saving...'])
    save(fullfile(saveLoc, [dNames{d} '_noMin_dat_S' num2str(sess)]),...
        'sl_rdm', 'rl_rdm', 'par','w', '-v7.3') % 'tsl', 'trl'
    
end

load(fullfile(saveLoc, 'sig_allSessions'))


clear newW
doNorm = 1;
spkSL = []; spkRL = [];
for sess = 1 : 8
    
    disp(['Session ' num2str(sess)])
    %% Combine across sessions
    
    clear rl_map rl_rdm sl_map sl_rdm sl_view par w
    load(fullfile(saveLoc, [dNames{d} '_noMin_dat_S' num2str(sess)]))
    if size(w, 1) > size(w, 2); w = w'; end
    
    newW{sess} = w;
    sn = 1; 
    
    newSL = sn * sl_rdm;
    newRL = sn * rl_rdm;
    
    if doNorm == 1
        mm = [];
        for ch = 0 : 1
            c0 = find(par.beh.cho_trg == ch);
            add = conv(squeeze(nanmean(newSL(c0, :), 1)), fw, 'same');
            mm = [mm add(1, tsl)];
            add = conv(squeeze(nanmean(newRL(c0, :), 1)), fw, 'same');
            mm = [mm add(1, trl)];
        end
        maxM = max(mm);
        minM = min(mm);
        exc = maxM - minM;
        
        normDim_paraS{d}{sess}.mini = minM;
        normDim_paraS{d}{sess}.maxi = maxM;
        normDim_paraS{d}{sess}.excur = exc;
        normDim_paraS{d}{sess}.mm = mm;
        
        newSL = (newSL - minM) ./ exc;
        newRL = (newRL - minM) ./ exc;
        
        normDim_sl{d}{sess} = newSL;
        normDim_rl{d}{sess} = newRL;
    end
    
    
    
    if nansum(nansum(newSL)) == 0
        newSL = nan(size(newSL));
    end
    spkSL = cat(1, spkSL, newSL);
    
    
    if nansum(nansum(newRL)) == 0
        newRL = nan(size(newRL));
    end
    spkRL = cat(1, spkRL, newRL);
    
end
sigSL.ramp_noMin = spkSL;
sigRL.ramp_noMin = spkRL;
sigSL.w.ramp_noMin = newW;

%% Plot

plotSig = {'ramp_noMin'};
sigName = {'S^{ramp w/o Min} [a.u.]'};
dimms = [63 20 61];

% From plot_singleTrials_highlights 


% For highlights plot
sigCohs = nanunique(sigSL.sig_coh); plotCohs = [6 2];
% yLims = {[0.3 2], [.3 2.2], [5 75]}; % Ramp, PC1, TinC
yLims_BLsub = {[-.6 1.4], [-.6 1.4], [-25 40]}; % Ramp, PC1, TinC
yLims_raw = {[-.5 1.3], [-.5 1.3], [0 50]}; % Ramp, PC1, TinC
yTicks = {[-.5 0 .5 1], [-.5 0 .5 1], [-20 0 20 40]};  % Ramp, PC1, TinC

fWin = 50;
plotFiltST = gausswin(fWin, 1.5); plotFiltST = plotFiltST / sum(plotFiltST);
xlims{1} = [0.2 0.5];
% tBL = xlims{1}(1) + fWin/70/1000 * [-1 1]; tBLv = find(par_sess{1}.tsl >= tBL(1) & par_sess{1}.tsl < tBL(2) );
tBLv = find(round(par_sess{1}.tsl, 3) == round(xlims{1}(1), 3)); 

% Trials to highlight
uTRL{2} = [1 3 19 20 21 65 45];
uTRL{6} = [2 3 5 6 11 14 16 17 65];
% Trials not to display
exclTRL{2} = [17 79];
exclTRL{6} = [4 24 52];
pltDash = {'-', [], ':'};

numTr = 80;
choiceCol = 1; choiceCols = {[0 0 0], [201 178 4]/255};
cohText{2} = {'25.6% contra'};
cohText{6} = {'0%'};










