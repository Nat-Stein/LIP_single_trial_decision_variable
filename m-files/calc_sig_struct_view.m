% calc_sig_struct_view

global dNames 
% dNames = {'TinC', 'TinI', 'Ramp', 'PC1', 'WhatD', 'WhenD', 'MinC', 'MinI'};

disp(['Finished generating signals until ' dNames{d}])

par = par_sess{1};
fw = ones(1, 50)/50;
tsl = find(par.tsl >= 0 & par.tsl <=0.6);
    
sig_coh = [];
session = [];
sigSL = [];
sigRL = [];
for d = 1 : length(dNames)
    spkSL = []; 
    
    dim = dNames{d};
    
    disp(['Merging signal structure viewing: ' dNames{d} '...'])
    
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
        load(fullfile(saveLoc, [dNames{d} '_dat_S' num2str(sess)]), 'sl_view', 'par')
        clear w 
        load(fullfile(saveLoc, [dNames{d} '_w_S' num2str(sess)]), 'w')
        if size(w, 1) > size(w, 2); w = w'; end
        
        newW{sess} = w;
        sn = 1; % Plot all and determine whether any dimension needs to be inverted!
        
        newSL = sn * sl_view;
        
        if doNorm == 1
            mm = [];
            for ch = [-1  1]
                c0 = find(sign(par.behView.sig_coh) == ch);
                add = conv(squeeze(nanmean(newSL(c0, :), 1)), fw, 'same');
                mm = [mm add(1, tsl)];
            end
            maxM = max(mm);
            minM = min(mm);
            exc = maxM - minM;
            
            normDim_paraS_view{d}{sess}.mini = minM;
            normDim_paraS_view{d}{sess}.maxi = maxM;
            normDim_paraS_view{d}{sess}.excur = exc;
            normDim_paraS_view{d}{sess}.mm = mm;
            
            newSL = (newSL - minM) ./ exc;
        end
        
        
        
        if nansum(nansum(newSL)) == 0
            newSL = nan(size(newSL));
        end
        spkSL = cat(1, spkSL, newSL);
        
        if d == 1
            nc = par.behView.sig_coh';
            sig_coh = cat(1, sig_coh, nc);
            n = sess * ones(size(nc));
            session = cat(1, session, n);
        end
        
    end
    switch dim
        case 'TinC'
            sigSL.TinC = spkSL;
            sigSL.w.TinC = newW;
            disp(num2str(dim))
        case 'TinI'
            sigSL.TinI = spkSL;
            sigSL.w.TinI = newW;
            disp(num2str(dim))
        case 'WhenD'
            sigSL.whenC = spkSL;
            sigSL.w.whenC = newW;
            disp(num2str(dim))
        case 'PC1'
            sigSL.PC1 = spkSL;
            sigSL.w.PC1 = newW;
            disp(num2str(dim))
        case 'MinC'
            sigSL.MinC = spkSL;
            sigSL.w.MinC = newW;
            disp(num2str(dim))
        case 'MinI'
            sigSL.MinI = spkSL;
            sigSL.w.MinI = newW;
            disp(num2str(dim))
        case 'WhatD'
            sigSL.whatD = spkSL;
            sigSL.w.whatD = newW;
            disp(num2str(dim))
        case 'Ramp_noLasso'
            sigSL.ramp_noLasso = spkSL;
            sigSL.w.ramp_noLasso = newW;
            disp(num2str(dim))
        case 'Ramp'
            sigSL.ramp = spkSL;
            sigSL.w.ramp = newW;
            disp(num2str(dim))
        case 'WhenD_Tin'
            sigSL.whenC_tin = spkSL;
            sigSL.w.whenC_tin = newW;
            disp(num2str(dim))
        case 'WhenD_noTin'
            sigSL.whenC_noTin = spkSL;
            sigSL.w.whenC_noTin = newW;
            disp(num2str(dim))
    end
end
sigSL.sig_coh = sig_coh;
sigSL.session = session;
sigSL.t = par.tsl; 

save(fullfile(saveLoc, 'sig_view_allSessions'),...
            'sigSL', 'tsl', 'fw', '-v7.3') 
