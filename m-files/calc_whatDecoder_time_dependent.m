% LIP choice decoder including/excluding Tin, Tout, Min populations -
% time-dependent <--
% Renames from calc_whatDecoder_acrossTPTs

load(fullfile(saveLoc, 'par_sess'))

winSize = 0.05 ; 
winCenter = 0.001 : 0.005 : 0.6;
winStart = winCenter - winSize/2;
winEnd = winCenter + winSize/2;
winCenter = mean([winStart ; winEnd]);

% Compute decoders including and excluding neural populations
exclName = {'All', 'noTin', 'noTout', 'noMinC', 'noTiTo', 'noTiMiC', 'noToMiC', 'noTiToMiC',...
    'justTin', 'justTout', 'justMinC', 'justTiTo', 'justTinMin', 'justMinCI', 'noTinMin', 'noMinCI'};


i_proj = find(round(winCenter, 3) == 0.451);

clear pc_dec se_dec w_dec

for sess = 1 : 8
    
    disp(['Session ' num2str(sess) ': loading normalized data...'])
    
    %% Load normalized activity for this session
    clear spkZnan par
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZnan', 'par')
    
    par = par_sess{sess};
    %% Set up decoder settings
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt > 0.6; % l = rt > 0.6  & coh >= 0;
    nL = sum(l);
    yChoice = par.beh.cho_trg(l);
    
    allNs = 1 : par.numCells;
    a = ismember(allNs, par.unit_Tin);
    notTin = allNs(a == 0);
    a = ismember(allNs, par.unit_Tout);
    notTout = allNs(a == 0);
    a = ismember(allNs, par.unit_MinC);
    notMinC = allNs(a == 0);
    a = ismember(allNs, [par.unit_Tin par.unit_Tout]);
    notTinTout = allNs(a == 0);
    a = ismember(allNs, [par.unit_Tin par.unit_Tout par.unit_MinC par.unit_MinI]);
    notTinMin = allNs(a == 0);
    a = ismember(allNs, [par.unit_MinC par.unit_MinI]);
    notMinCI = allNs(a == 0);
    
    excl = {[],...                          % all neurons (Fig 4a)
        par.unit_Tin,...                    % all neurons but TinC
        par.unit_Tout,...                   % all neurons but TinI
        par.unit_MinC,...                   % all neurons but MinC
        [par.unit_Tin par.unit_Tout],...    % all neurons but TinC and TinI (Fig S9)
        [par.unit_Tin par.unit_MinC],...    % all neurons but TinC and MinC
        [par.unit_Tout par.unit_MinC],...   % all neurons but TinI and MinC
        [par.unit_Tin par.unit_Tout par.unit_MinC],...  % all neurons but TinC, TinI & MinC
        notTin,...                          % only TinC      
        notTout,...                         % only TinI
        notMinC,...                         % only MinC
        notTinTout,...                      % only TinC & TinI  (Fig S9)
        notTinMin,...                       % only TinC & MinC
        notMinCI,...                        % only MinC & MinI
        [par.unit_Tin par.unit_Tout par.unit_MinC par.unit_MinI],... % all but TinC, TinI, MinC & MinI  (Fig S9)
        [par.unit_MinC par.unit_MinI],...   % all but MinC & MinI
        };

    %% Compute what decoders
    disp(['Session ' num2str(sess) ': Compute What decoders...'])
    
    excl_sess{sess} = excl;
    
    for e = 1 : length(excl) % loop through exclusion criteria
        
        disp(['Condition ' num2str(e) '/' num2str(length(excl))])
        
        clear pcDD_test pDDtest_se useN
        
        useN = ones(par_sess{sess}.numCells, 1); useN(excl{e}) = 0; useN(find(spkZnan.s==0)) = 0;
        useN_dec{e, sess} = find(useN==1); % variable useN_dec stores the indices of neurons included in each decoder
        
        if sum(useN) > 0
            for i = 1 : length(winStart)
                
                dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
                
                xDD = squeeze(nanmean(spkZnan.SL(useN==1, l, dW), 3));
                
                xDD = xDD - nanmean(xDD,2) ;
                lDD = (nanmean(~isnan(xDD), 1) == 1)';
                
                train = ~logical(mod(1:nL,2))' ;
                xD = xDD(:, train&lDD)' ;
                xTest = xDD(:, ~train&lDD)' ;
                c = yChoice(train&lDD);
                
                % Fit logistic
                [wDD,sDD]     =  lassoglm(-xD, c', 'binomial','link','logit','Lambda',0.03);
                wDD = [0 ; wDD] ;
                
                choicePredTrain  = glmval(wDD,  -xD, 'logit') > 0.5;
                
                pcDD_train  =  sum(choicePredTrain'==yChoice(train&lDD)) ./ sum(train&lDD) ;
                
                choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
                
                pcDD_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                pDDtest_se(i) = sqrt( (pcDD_test(i)*(1-pcDD_test(i))) / sum(~train&lDD)) ;
                
                w_dec{e, sess}(i, :) = wDD;
                
            end
            pc_dec{e, sess} = pcDD_test;
            se_dec{e, sess} = pDDtest_se;
            
            if e == 1
                allNs = par_sess{sess}.numCells;
                w = ones(1, length(allNs));
                addWeights = w_dec{e, sess}(i_proj, 2:end);
                w(1, useN_dec{e, sess}) = addWeights; 
                ff = fullfile(saveLoc,...
                    [dNames{5} '_w_S' num2str(sess)]);
                save(ff, 'w')
            end
            
        end
    end
end

disp('Saving stimulus-locked decoders...')
save(fullfile(saveLoc, 'decoder_excl_SL'), 'w_dec', 'pc_dec',...
    'se_dec',...
    'excl_sess', 'exclName', 'useN_dec', 'winCenter', 'winStart', 'winEnd', 'winSize')

%% Compute time-dependent decoders on response-locked data

disp('Computing response-locked decoders...')

i_proj = find(round(winCenter, 3) == 0.451);

trl_dec = -0.4 : 0.005 : 0; wBin = 0.005;
calcWhenD = 0;

clear w_decR
for sess = 1 : 8
    
    disp(['Session ' num2str(sess) ': loading normalized response-locked data...'])
    clear spkZnan par
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZnan')
    par = par_sess{sess};
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt > 0.6; % l = rt > 0.6  & coh >= 0;
    nL = sum(l);
    yChoice = par.beh.cho_trg(l);
    
    excl = excl_sess{sess};
    
    disp(['Session ' num2str(sess) ': fitting decoders...'])
    
    for e = 1 : length(excl)
        
        disp([num2str(e) '/' num2str(length(excl))])
        
        clear pcDD_test pDDtest_se useN
        
        useN = useN_dec{e, sess};
        
        if sum(useN) > 0
            w_proj = w_dec{e, sess}(i_proj, :);
            
            for i = 1 : length(trl_dec)
                
                dW = find(par.trl >= trl_dec(i)-winSize/2 & par.trl <= trl_dec(i)+winSize/2); % Decoding window relative to dots onset
                
                xDD = squeeze(nanmean(spkZnan.RL(useN, l, dW), 3));
                
                xDD = xDD - nanmean(xDD,2) ;
                lDD = (nanmean(~isnan(xDD), 1) == 1)';
                
                train = ~logical(mod(1:nL,2))' ;
                xD = xDD(:, train&lDD)' ;
                xTest = xDD(:, ~train&lDD)' ;
                c = yChoice(train&lDD);
                
                [wDD,sDD]     =  lassoglm(-xD, c', 'binomial','link','logit','Lambda',0.03);
                wDD = [0 ; wDD] ;
                
                choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
                
                pcDD_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)') ./ sum(~train&lDD) ;
                
                pDDtest_se(i) = sqrt( (pcDD_test(i)*(1-pcDD_test(i))) / sum(~train&lDD)) ;
                
                w_decR{e, sess}(i, :) = wDD;
                
                choicePredTest  = glmval(w_proj',  -xTest, 'logit') > 0.5;
                
                pcDD450_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)') ./ sum(~train&lDD) ;
                
                pDDtest450_se(i) = sqrt( (pcDD450_test(i)*(1-pcDD450_test(i))) / sum(~train&lDD)) ;
                
            end
            pc_decR{e, sess} = pcDD_test;
            se_decR{e, sess} = pDDtest_se;
            
            pc450_decR{e, sess} = pcDD450_test;
            se450_decR{e, sess} = pDDtest450_se;
            
            if e == 1
                % Weights without Tin at same time point and weight without Tin from 450ms applied at each time point
                goodN = useN_dec{e, sess};
                isTin = zeros(1, size(spkZnan.RL, 1)); isTin(par.unit_Tin) = 1; isTin_n = isTin(goodN);
                useN = goodN(~isTin_n);
                isTiTo = zeros(1, size(spkZnan.RL, 1)); isTiTo([par.unit_Tin par.unit_Tout]) = 1; isTiTo_n = isTiTo(goodN);
                useN_tito = goodN(~isTiTo_n);
                
                % Weights without Tin at 450ms
                w_proj = w_dec{e, sess}(i_proj, 2:end);
                w_proj_nTin = w_proj(~isTin_n);
                w_proj_nTin = [0 w_proj_nTin];
                
                % Weights without Tin & Tout at 450ms
                w_proj = w_dec{e, sess}(i_proj, 2:end);
                w_proj_nTiTo = w_proj(~isTiTo_n);
                w_proj_nTiTo = [0 w_proj_nTiTo];
                
                for i = 1 : length(trl_dec)
                    
                    % Weights without Tin at same time point
                    w_proj = w_dec{e, sess}(i, 2:end);
                    w_proj_t = w_proj(~isTin_n);
                    w_proj_t = [0 w_proj_t];
                    
                    dW = find(par.trl >= trl_dec(i)-winSize/2 & par.trl <= trl_dec(i)+winSize/2);
                    dW_when = find(par.trl >= trl_dec(i)-winSize/2 & par.trl <= trl_dec(i)+winSize/2);
                    
                    % Without Tin
                    xDD = squeeze(nanmean(spkZnan.RL(useN, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    %xD = xDD(:, train&lDD)' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    % Weights without Tin at 450ms
                    choicePredTest  = glmval(w_proj_nTin',  -xTest, 'logit') > 0.5;
                    pcR_nTin_450_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)') ./ sum(~train&lDD) ;
                    pcR_nTin_450_test_se(i) = sqrt( (pcR_nTin_450_test(i)*(1-pcR_nTin_450_test(i))) / sum(~train&lDD)) ;
                    
                    % Weights without Tin at same time point
                    choicePredTest  = glmval(w_proj_t',  -xTest, 'logit') > 0.5;
                    pcR_nTin_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)') ./ sum(~train&lDD) ;
                    pcR_nTin_test_se(i) = sqrt( (pcR_nTin_test(i)*(1-pcR_nTin_test(i))) / sum(~train&lDD)) ;
                    
                    % ---------------------------------------------------
                    % Weights without Tin & Tout at same time point
                    w_proj = w_dec{e, sess}(i, 2:end);
                    w_proj_t = w_proj(~isTiTo_n);
                    w_proj_t = [0 w_proj_t];
                    
                    % Without Tin % Tout
                    xDD = squeeze(nanmean(spkZnan.RL(useN_tito, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    % Weights without Tin at 450ms
                    choicePredTest  = glmval(w_proj_nTiTo',  -xTest, 'logit') > 0.5;
                    pcR_nTiTo_450_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)') ./ sum(~train&lDD) ;
                    pcR_nTiTo_450_test_se(i) = sqrt( (pcR_nTiTo_450_test(i)*(1-pcR_nTiTo_450_test(i))) / sum(~train&lDD)) ;
                    
                    % Weights without Tin at same time point
                    choicePredTest  = glmval(w_proj_t',  -xTest, 'logit') > 0.5;
                    pcR_nTiTo_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)') ./ sum(~train&lDD) ;
                    pcR_nTiTo_test_se(i) = sqrt( (pcR_nTiTo_test(i)*(1-pcR_nTiTo_test(i))) / sum(~train&lDD)) ;
                    
                    
                    % ---------------------------------------------------
                    % With Tin
                    xDD = squeeze(nanmean(spkZnan.RL(goodN, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    %xD = xDD(:, train&lDD)' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    if calcWhenD == 1
                        % Weights of when-decoder using dimension traces
                        xDD_dim = nanmean(allDim_rl{7}{1, sess}(l, dW_when), 2); xDD_dim = xDD_dim - nanmean(xDD_dim);
                        tTrain_dim = xDD_dim(train&lDD)' ; tTrain_dim = [zeros(size(tTrain_dim)); tTrain_dim];
                        [wDD,sDD]     =  lassoglm(-tTrain_dim', c, 'binomial','link','logit','Lambda',0);
                        tTest_dim = xDD_dim(~train&lDD)' ;
                        choicePredTest  = glmval(wDD,  -tTest_dim, 'logit') > 0.5; % When decoder weights do not get applied to Z but 0 data!!!
                        pcR_whenC_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                        pcR_whenC_test_se(i) = sqrt( (pcR_whenC_test(i)*(1-pcR_whenC_test(i))) / sum(~train&lDD)) ;
                        %                     choicePredTest  = glmval(w_proj_whenC',  -xTest, 'logit') > 0.5;
                        %                     pcR_whenC_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                        %                     pcR_whenC_test_se(i) = sqrt( (pcR_whenC_test(i)*(1-pcR_whenC_test(i))) / sum(~train&lDD)) ;
                        
                        % Weights of when-decoder with just Tin neuronsxDD_dim = nanmean(allDim_rl{7}{1, sess}(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim);
                        xDD_dim = nanmean(allDim_rl{70}{1, sess}(l, dW_when), 2); xDD_dim = xDD_dim - nanmean(xDD_dim);
                        tTrain_dim = xDD_dim(train&lDD)' ; tTrain_dim = [zeros(size(tTrain_dim)); tTrain_dim];
                        [wDD,sDD]     =  lassoglm(-tTrain_dim', c, 'binomial','link','logit','Lambda',0);
                        tTest_dim = xDD_dim(~train&lDD)' ;
                        choicePredTest  = glmval(wDD,  -tTest_dim, 'logit') > 0.5; % When decoder weights do not get applied to Z but 0 data!!!
                        pcR_whenC_Tin_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                        pcR_whenC_Tin_test_se(i) = sqrt( (pcR_whenC_Tin_test(i)*(1-pcR_whenC_Tin_test(i))) / sum(~train&lDD)) ;
                        %                     choicePredTest  = glmval(w_proj_whenC_Tin',  -xTest, 'logit') > 0.5;
                        %                     pcR_whenC_Tin_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                        %                     pcR_whenC_Tin_test_se(i) = sqrt( (pcR_whenC_Tin_test(i)*(1-pcR_whenC_Tin_test(i))) / sum(~train&lDD)) ;
                        
                        % Weights of when-decoder without Tin neurons
                        xDD_dim = nanmean(allDim_rl{71}{1, sess}(l, dW_when), 2); xDD_dim = xDD_dim - nanmean(xDD_dim);
                        tTrain_dim = xDD_dim(train&lDD)' ; tTrain_dim = [zeros(size(tTrain_dim)); tTrain_dim];
                        [wDD,sDD]     =  lassoglm(-tTrain_dim', c, 'binomial','link','logit','Lambda',0);
                        tTest_dim = xDD_dim(~train&lDD)' ;
                        choicePredTest  = glmval(wDD,  -tTest_dim, 'logit') > 0.5; % When decoder weights do not get applied to Z but 0 data!!!
                        pcR_whenC_nTin_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                        pcR_whenC_nTin_test_se(i) = sqrt( (pcR_whenC_nTin_test(i)*(1-pcR_whenC_nTin_test(i))) / sum(~train&lDD)) ;
                        %                     choicePredTest  = glmval(w_proj_whenC_nTin',  -xTest, 'logit') > 0.5;
                        %                     pcR_whenC_nTin_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                        %                     pcR_whenC_nTin_test_se(i) = sqrt( (pcR_whenC_nTin_test(i)*(1-pcR_whenC_nTin_test(i))) / sum(~train&lDD)) ;
                    end
                end
                pc450_Tin0_decR{e, sess} = pcR_nTin_450_test;
                se450_Tin0_decR{e, sess} = pcR_nTin_450_test_se;
                
                pc_Tin0_decR{e, sess} = pcR_nTin_test;
                se_Tin0_decR{e, sess} = pcR_nTin_test_se;
                
                pc450_TiTo0_decR{e, sess} = pcR_nTiTo_450_test;
                se450_TiTo0_decR{e, sess} = pcR_nTiTo_450_test_se;
                
                pc_TiTo0_decR{e, sess} = pcR_nTiTo_test;
                se_TiTo0_decR{e, sess} = pcR_nTiTo_test_se;
                
                if calcWhenD == 1
                    pc_whenC_decR{e, sess} = pcR_whenC_test;
                    se_whenC_decR{e, sess} = pcR_whenC_test_se;
                    
                    pc_whenC_Tin_decR{e, sess} = pcR_whenC_Tin_test;
                    se_whenC_Tin_decR{e, sess} = pcR_whenC_Tin_test_se;
                    
                    pc_whenC_nTin_decR{e, sess} = pcR_whenC_nTin_test;
                    se_whenC_nTin_decR{e, sess} = pcR_whenC_nTin_test_se;
                end
            end
        end
    end
end

save(fullfile(saveLoc, 'decoder_excl_RL'), 'w_decR', 'pc_decR',...
    'se_decR', 'excl_sess', 'exclName', 'useN_dec', 'trl_dec', 'winCenter', 'winStart', 'winEnd', 'winSize')

clear pc_whenC_sign_decR se_whenC_sign_decR
for sess = 1 : 8
    
    load(fullfile(saveLoc, [dNames{d} '_dat_S' num2str(sess)]),...
        'sl_rdm', 'rl_rdm', 'par')
    par = par_sess{sess};
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt > 0.6;
    nL = sum(l);
    yChoice = par.beh.cho_trg(l);
    
    clear pcR_whenC_sign_test seR_whenC_sign_test
    for i = 1 : length(trl_dec)
        
        dW = find(par.trl >= trl_dec(i)-winSize/2 & par.trl <= trl_dec(i)+winSize/2);
        % When-decoder sign test
        xDD_dim = nanmean(rl_rdm(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); lDD = (~isnan(xDD_dim) == 1)';
        cho_sign = sign(xDD_dim);
        choicePredSign = cho_sign(lDD) == -1;
        pcR_whenC_sign_test(i)  =  sum(choicePredSign==yChoice(lDD)') ./ sum(lDD) ;
        
        seR_whenC_sign_test(i) = sqrt( (pcR_whenC_sign_test(i)*(1-pcR_whenC_sign_test(i))) / sum(lDD)) ;
        
    end
    pc_whenC_sign_decR{1, sess} = pcR_whenC_sign_test;
    se_whenC_sign_decR{1, sess} = seR_whenC_sign_test;
    
    
end

save(fullfile(saveLoc, 'decoder_excl_RL'), 'w_dec', 'pc_dec', 'se_dec',...
    'excl_sess', 'excl', 'exclName', 'useN_dec', 'winSize', 'trl_dec',...
    'pc_decR', 'se_decR', 'pc450_decR', 'se450_decR',...
    'pc450_Tin0_decR', 'se450_Tin0_decR', 'pc_Tin0_decR', 'se_Tin0_decR',...
    'pc450_TiTo0_decR', 'se450_TiTo0_decR', 'pc_TiTo0_decR', 'se_TiTo0_decR',...
    'pc_whenC_sign_decR', 'se_whenC_sign_decR')
