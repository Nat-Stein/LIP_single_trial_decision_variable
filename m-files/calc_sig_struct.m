% calc_sig_struct

global dNames 
% dNames = {'TinC', 'TinI', 'Ramp', 'PC1', 'WhatD', 'WhenD', 'MinC', 'MinI'};

for d = 1 : length(dNames)
    
    dim = dNames{d};
    
    disp(['Generating signals: ' dNames{d} '...'])
    
    switch dim
        case 'PC1'
            nDat = 1;
        case 'WhatD'
            nDat = 1;
        case 'Ramp'
            nDat = 2;
        otherwise
            nDat = 0;
    end
    
    for sess = 1 : 8
        
        disp([dNames{d} ', Session ' num2str(sess)])
        
        clear w 
        ff = fullfile(saveLoc,...
            [dNames{d} '_w_S' num2str(sess)]);
        load(ff, 'w')
        if size(w, 1) > size(w, 2); w = w'; end
        
        clear par dd beh behMap behView spk spkMap spkView media stdev
        clear spkZ spkMapZ spkViewZ spkZnan spkMapZnan spkViewZnan
        cLen = 200;
        if nDat == 0 
            %% Use non-normalized data  ---------------------------------
            disp(['Session ' num2str(sess) ': loading...'])
            load(fullfile(saveLoc, ['spk_S' num2str(sess)]))
            
            disp(['Session ' num2str(sess) ': projecting data...'])
            tpts = 1 : size(spk.SLnan, 3);
            tpts_r = 1 : size(spk.RLnan, 3);
            % RDM: Stimulus-locked data
            Cin = squeeze(spk.SLnan(:, :, tpts)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            sl_rdm = C;
            % RDM: Response-locked data
            Cin = squeeze(spk.RLnan(:, :, tpts_r)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            rl_rdm = C;
            % Mapping: Stimulus-locked data
            Cin = squeeze(spkMap.SLnan(:, :, tpts)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            sl_map = C;
            % Mapping: Response-locked data
            Cin = squeeze(spkMap.RLnan(:, :, tpts_r)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            rl_map = C;
            % Viewing: Stimulus-locked data
            Cin = squeeze(spkView.SLnan(:, :, tpts)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            sl_view = C;
            clear Cin C
            
        elseif  nDat == 1 
            %% Use normalized data ---------------------------------
            % PC1 and What-decoder
            disp(['Session ' num2str(sess) ': loading...'])
            load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]), 'spkZnan', 'spkMapZnan', 'spkViewZnan', 'par')
            
            disp(['Session ' num2str(sess) ': projecting data...'])
            w_orig = w;
            if length(w) > length(spkZnan.goodN); w = w(spkZnan.goodN); end
            tpts = 1 : size(spkZnan.SL, 3);
            tpts_r = 1 : size(spkZnan.RL, 3);
            % RDM: Stimulus-locked data
            Cin = squeeze(spkZnan.SL(spkZnan.goodN, :, tpts)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            sl_rdm = C;
            % RDM: Response-locked data
            Cin = squeeze(spkZnan.RL(spkZnan.goodN, :, tpts_r)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            rl_rdm = C;
            % Mapping: Stimulus-locked data
            Cin = squeeze(spkMapZnan.SL(spkMapZnan.goodN, :, tpts)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            sl_map = C;
            % Mapping: Response-locked data
            Cin = squeeze(spkMapZnan.RL(spkMapZnan.goodN, :, tpts_r)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            rl_map = C;
            % Viewing: Stimulus-locked data
            Cin = squeeze(spkViewZnan.SL(spkViewZnan.goodN, :, tpts)); C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w, Cin(:, :, tps)));
            end
            sl_view = C;
            w = w_orig;
            
        else
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
            sl_rdm = C; figure; imagesc(C)
            % RDM: Response-locked data
            Cin = (squeeze(spk.RLnan(goodN, :, :)) - mMat) ./ sMat;
            C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w(goodN), Cin(:, :, tps)));
            end
            rl_rdm = C;
            
            % Mapping: Stimulus-locked data
            mMat = repmat(media(goodN)', 1, size(spkMap.SLnan, 2), size(spkMap.SLnan, 3));
            sMat = repmat(stdev(goodN)', 1, size(spkMap.SLnan, 2), size(spkMap.SLnan, 3));
            Cin = (squeeze(spkMap.SLnan(goodN, :, :)) - mMat) ./ sMat;
            C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w(goodN), Cin(:, :, tps)));
            end
            sl_map = C;
            % Mapping: Response-locked data
            Cin = (squeeze(spkMap.RLnan(goodN, :, :)) - mMat) ./ sMat; C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w(goodN), Cin(:, :, tps)));
            end
            rl_map = C;
            % Viewing: Stimulus-locked data
            mMat = repmat(media(goodN)', 1, size(spkView.SLnan, 2), size(spkView.SLnan, 3));
            sMat = repmat(stdev(goodN)', 1, size(spkView.SLnan, 2), size(spkView.SLnan, 3));
            Cin = (squeeze(spkView.SLnan(goodN, :, :)) - mMat) ./ sMat;
            C = nan(size(Cin, 2), size(Cin, 3));
            for cc = 1 : ceil(size(Cin, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(Cin, 3))];
                C(:, tps) = squeeze(mmx('mult', w(goodN), Cin(:, :, tps)));
            end
            sl_view = C;
        end
        
        %% Save single session data matrices
        disp(['Session ' num2str(sess) ': saving...'])
        save(fullfile(saveLoc, [dNames{d} '_dat_S' num2str(sess)]),...
            'sl_rdm', 'rl_rdm', 'sl_map', 'rl_map', 'sl_view', 'par','w', '-v7.3') % 'tsl', 'trl'
        clear sl_map rl_map sl_view sl_rdm rl_rdm 
        clear beh behMap behView spk spkMap spkView
        clear spkZ spkMapZ spkViewZ spkZnan spkMapZnan spkViewZnan
        disp(['Session ' num2str(sess) ': completed.'])
    end
end

disp(['Finished generating signals until ' dNames{d}])

fw = ones(1, 50)/50;
tsl = find(par.tsl >= 0 & par.tsl <=0.6);
trl = find(par.trl >= -0.6 & par.trl <=0);
    
sig_coh = [];
correct = [];
session = [];
rt = [];
choice = [];
sigSL = [];
sigRL = [];
for d = 1 : length(dNames)
    spkSL = []; spkRL = [];
    
    dim = dNames{d};
    
    disp(['Merging signal structure: ' dNames{d} '...'])
    
    switch dim
        case 'PC1'
            doNorm = 1;
        case 'WhatD'
            doNorm = 1;
        case 'Ramp'
            doNorm = 1;
        case 'WhenD'
            doNorm = 1;
        otherwise
            doNorm = 0;
    end
    
    clear newW
    for sess = 1 : 8
        
        disp(['Session ' num2str(sess)])
        %% Combine across sessions
        
        clear rl_map rl_rdm sl_map sl_rdm sl_view par
        load(fullfile(saveLoc, [dNames{d} '_dat_S' num2str(sess)]))
        clear w 
        load(fullfile(saveLoc, [dNames{d} '_w_S' num2str(sess)]), 'w')
        if size(w, 1) > size(w, 2); w = w'; end
        
        newW{sess} = w;
        sn = 1; % Plot all and determine whether any dimension needs to be inverted!
        
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
        
        if d == 1
            nc = par.beh.sig_coh';
            sig_coh = cat(1, sig_coh, nc);
            nc = par.beh.correct';
            correct = cat(1, correct, nc);
            nc = par.beh.rt';
            rt = cat(1, rt, nc);
            nc = par.beh.cho_trg';
            choice = cat(1, choice, nc);
            n = sess * ones(size(nc));
            session = cat(1, session, n);
        end
        
    end
    switch dim
        case 'TinC'
            sigSL.TinC = spkSL;
            sigRL.TinC = spkRL;
            sigSL.w.TinC = newW;
            disp(num2str(dim))
        case 'TinI'
            sigSL.TinI = spkSL;
            sigRL.TinI = spkRL;
            sigSL.w.TinI = newW;
            disp(num2str(dim))
        case 'WhenD'
            sigSL.whenC = spkSL;
            sigRL.whenC = spkRL;
            sigSL.w.whenC = newW;
            disp(num2str(dim))
        case 'PC1'
            sigSL.PC1 = spkSL;
            sigRL.PC1 = spkRL;
            sigSL.w.PC1 = newW;
            disp(num2str(dim))
        case 'MinC'
            sigSL.MinC = spkSL;
            sigRL.MinC = spkRL;
            sigSL.w.MinC = newW;
            disp(num2str(dim))
        case 'MinI'
            sigSL.MinI = spkSL;
            sigRL.MinI = spkRL;
            sigSL.w.MinI = newW;
            disp(num2str(dim))
        case 'WhatD'
            sigSL.whatD = spkSL;
            sigRL.whatD = spkRL;
            sigSL.w.whatD = newW;
            disp(num2str(dim))
        case 'Ramp_noLasso'
            sigSL.ramp_noLasso = spkSL;
            sigRL.ramp_noLasso = spkRL;
            sigSL.w.ramp_noLasso = newW;
            disp(num2str(dim))
        case 'Ramp'
            sigSL.ramp = spkSL;
            sigRL.ramp = spkRL;
            sigSL.w.ramp = newW;
            disp(num2str(dim))
        case 'WhenD_Tin'
            sigSL.whenC_tin = spkSL;
            sigRL.whenC_tin = spkRL;
            sigSL.w.whenC_tin = newW;
            disp(num2str(dim))
        case 'WhenD_noTin'
            sigSL.whenC_noTin = spkSL;
            sigRL.whenC_noTin = spkRL;
            sigSL.w.whenC_noTin = newW;
            disp(num2str(dim))
    end
end
sigSL.sig_coh = sig_coh;
sigSL.correct = correct;
sigSL.session = session;
sigSL.rt = rt;
sigSL.choice = choice;
sigSL.t = par.tsl; 

sigRL.sig_coh = sig_coh;
sigRL.correct = correct;
sigRL.session = session;
sigRL.rt = rt;
sigRL.choice = choice;
sigRL.t = par.trl; 
sigRL.w = sigSL.w;

save(fullfile(saveLoc, 'sig_allSessions'),...
            'sigSL', 'sigRL', 'tsl', 'trl', 'fw', '-v7.3') 
