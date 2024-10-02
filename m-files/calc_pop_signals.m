% calc_pop_signals


% load(fullfile(saveLoc, 'par_sess'))

for sess = 1 : 8
    
    clear par dd beh behMap behView spk spkMap spkView
    clear spkZ spkMapZ spkViewZ spkZnan spkMapZnan spkViewZnan
    
    load(fullfile(saveLoc, ['spk_S' num2str(sess)]))
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]))
    
    for d = 1 : length(dNames)
        clear w 
        dim = dNames{d};
        ff = fullfile(saveLoc,...
            [dNames{d} '_w_S' num2str(sess)]);
        load(ff, 'w')
        if size(w, 1) > size(w, 2); w = w'; save(ff, 'w'); end
        
        % Does the dimension use normalized data?
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
        % Use z-scored data for ramp direction!
        
        clear sl_rdm rl_rdm sl_map rl_map sl_view
        if nDat == 1 % Use normalized data
            % RDM: Stimulus-locked data
            Cin = squeeze(spkZnan.SL(spkZnan.goodN, :, tpts));
            sl_rdm = squeeze(mmx('mult', w(spkZnan.goodN), Cin));
            % RDM: Response-locked data
            Cin = squeeze(spkZnan.RL(spkZnan.goodN, :, tpts_r));
            rl_rdm = squeeze(mmx('mult', w(spkZnan.goodN), Cin));
            % Mapping: Stimulus-locked data
            Cin = squeeze(spkMapZnan.SL(spkMapZnan.goodN, :, tpts));
            sl_map = squeeze(mmx('mult', w(spkMapZnan.goodN), Cin));
            % Mapping: Response-locked data
            Cin = squeeze(spkMapZnan.RL(spkMapZnan.goodN, :, tpts_r));
            rl_map = squeeze(mmx('mult', w(spkMapZnan.goodN), Cin));
            % Viewing: Stimulus-locked data
            Cin = squeeze(spkViewZnan.SL(spkViewZnan.goodN, :, tpts));
            sl_view = squeeze(mmx('mult', w(spkViewZnan.goodN), Cin));
        elseif nDat == 2 % Use z-scored data (ramp)
            
            % RDM: Stimulus-locked data
            Cin = squeeze(spkZnan.SL(spkZnan.goodN, :, tpts));
            sl_rdm = squeeze(mmx('mult', w(spkZnan.goodN), Cin));
            % RDM: Response-locked data
            Cin = squeeze(spkZnan.RL(spkZnan.goodN, :, tpts_r));
            rl_rdm = squeeze(mmx('mult', w(spkZnan.goodN), Cin));
            % Mapping: Stimulus-locked data
            Cin = squeeze(spkMapZnan.SL(spkMapZnan.goodN, :, tpts));
            sl_map = squeeze(mmx('mult', w(spkMapZnan.goodN), Cin));
            % Mapping: Response-locked data
            Cin = squeeze(spkMapZnan.RL(spkMapZnan.goodN, :, tpts_r));
            rl_map = squeeze(mmx('mult', w(spkMapZnan.goodN), Cin));
            % Viewing: Stimulus-locked data
            Cin = squeeze(spkViewZnan.SL(spkViewZnan.goodN, :, tpts));
            sl_view = squeeze(mmx('mult', w(spkViewZnan.goodN), Cin));
            
            
        else % Use non-normalized data
            % RDM: Stimulus-locked data
            Cin = squeeze(spk.SLnan(:, :, tpts));
            sl_rdm = squeeze(mmx('mult', w, Cin));
            % RDM: Response-locked data
            Cin = squeeze(spk.RLnan(:, :, tpts_r));
            rl_rdm = squeeze(mmx('mult', w, Cin));
            % Mapping: Stimulus-locked data
            Cin = squeeze(spkMap.SLnan(:, :, tpts));
            sl_map = squeeze(mmx('mult', w, Cin));
            % Mapping: Response-locked data
            Cin = squeeze(spkMap.RLnan(:, :, tpts_r));
            rl_map = squeeze(mmx('mult', w, Cin));
            % Viewing: Stimulus-locked data
            Cin = squeeze(spkView.SLnan(:, :, tpts));
            sl_view = squeeze(mmx('mult', w, Cin));
        end
        
        save(fullfile(saveLoc, [dNames{d} '_dat_S' num2str(sess)]),...
            'sl_rdm', 'rl_rdm', 'sl_map', 'rl_map', 'sl_view', 'tsl', 'trl', '-v7.3')
        
    end
end