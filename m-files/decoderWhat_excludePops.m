



%% Grab weights at 450ms and make them into dimensions - stimulus aligned
% Add dimension with same weights but Tin set to zero
% I'm not using this, right???

i_proj = find(round(winCenter, 3) == 0.451);

for sess = 1 : 8
    par = par_sess{sess};
    
    % Load normalized activity for this session
    clear spkZ
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZ')
    
    clear useDat
    for tp = 1:length(winStart)
        tpts = find(par.tsl >= winStart(tp) & par.tsl <= winEnd(tp));
        useDat(:, :, tp) = nanmean(spkZ.SL(:, :, tpts), 3);
    end
    
    useN = ones(length(allNs), 1); useN(excl{e}) = 0; useN(spkZ.s==0) = 0;
    useN_dec{e, sess} = find(useN==1);
    
    for e = 1 : size(pc_dec, 1)
        disp(['Sess ' num2str(sess) 'e ' num2str(e) '/' num2str(length(excl))])
        
        if size(w_dec{e, sess}, 1) >0
            w = w_dec{e, sess}(i_proj, 2:end); % used to be: w = w_dec{e, sess}(i_proj, :);
            %goodN = ones(size(rl_sessMatZ_short{sess}, 1), 1); goodN(excl{e}) = 0; goodN(sessZ{sess}.s==0) = 0;
            goodN = useN_dec{e, sess};
            
            C = squeeze(mmx('mult', w, useDat(goodN, :, :)));
            
            sl_exDec{e} = C;
            
            if e == 1
                % Set Tin neuron's activity to zero
                w = w_dec{e, sess}(i_proj, 2:end);
                goodN = useN_dec{e, sess};
                isTin = zeros(1, size(useDat, 1)); isTin(par.unit_Tin) = 1; isTin_n = isTin(goodN);
                dat = useDat(goodN, :, :);
                
                C = squeeze(mmx('mult', w(~isTin_n), dat(~isTin_n, :, :)));
                
                sl_exDec_Tin0{e} = C;
            end
            
        end
    end
    save(fullfile(saveLoc, ['sl_decoder' num2str(sess)]),...
        'sl_exDec', 'sl_exDec_Tin0', 'winStart', 'winCenter', 'winEnd', 'se_dec', 'excl', 'exclName', 'i_proj')
end

% save('E:\Matalb analyses\sl_decoder_excluded_230522', 'sl_exDec', 'sl_exDec_Tin0', 'se_dec', 'excl', 'exclName')

%% Decoder performace: stimulus locked
%% 1. weights from 450ms applied to all time points
%% 2. weights from 450ms with Tin neurons set to zero applied to all time points
%% 3. weights with Tin neurons set to zero applied to same point
%% 4. weights of when decoder applied to all time points
%% 5. weights of when decoder with just Tin neurons applied to all time points
%% 6. weights of when decoder without Tin neurons applied to all time points

% Need to load weights of other dimensions?
load('E:\Matalb analyses\sl_allDim_minimal_230522', 'allW', 'allDim_sl', 'allDim_rl')
load(fullfile(par_sess{1}.saveLoc, 'decoder_excl_SL'))
i_proj = find(round(winCenter, 3) == 0.451);

% Which of these outputs do I actually need????
for sess = 1 : 8
    par = par_sess{sess};
    
    % Load normalized activity for this session
    clear spkZ
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZ')
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt>0.6; % l = rt>0.6  & coh >= 0;
    nL = sum(l);
    yChoice = par.beh.cho_trg(l);
    
    for e = 1 : length(excl_sess{sess})
        
        disp(['Sess ' num2str(sess) 'e ' num2str(e) '/' num2str(length(excl))])
        
        clear pcDD_test pDDtest_se useN pDDtest550_se pcDD550_test
        
        useNs = useN_dec{e, sess};
        
        if size(w_dec{e, sess}, 1) > 0
            
            w_proj = w_dec{e, sess}(i_proj, :);
            for i = 1 : length(winStart)
                
                dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
                
                
                xDD = squeeze(nanmean(spkZ.SL(useNs, l, dW), 3));
                
                xDD = xDD - nanmean(xDD,2) ;
                lDD = (nanmean(~isnan(xDD), 1) == 1)';
                
                train = ~logical(mod(1:nL,2))' ;
                %xD = xDD(:, train&lDD)' ;
                xTest = xDD(:, ~train&lDD)' ;
                c = yChoice(train&lDD);
                
                choicePredTest  = glmval(w_proj',  -xTest, 'logit') > 0.5;
                
                pcDD450_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                pDDtest450_se(i) = sqrt( (pcDD450_test(i)*(1-pcDD450_test(i))) / sum(~train&lDD)) ;
                
            end
            pc450_dec{e, sess} = pcDD450_test;
            se450_dec{e, sess} = pDDtest450_se;
            
            
            if e == 1 
                % Weights without Tin at same time point and weight without Tin from 450ms applied at each time point
                goodN = useN_dec{e, sess};
                isTin = zeros(1, size(spkZ.SL, 1)); isTin(par.unit_Tin) = 1; isTin_n = isTin(goodN);
                useN = goodN(~isTin_n);
                isTiTo = zeros(1, size(spkZ.SL, 1)); isTiTo([par.unit_Tin par.unit_Tout]) = 1; isTiTo_n = isTiTo(goodN);
                useN_tito = goodN(~isTiTo_n);
                
                % Weights without Tin at 450ms
                w_proj = w_dec{e, sess}(i_proj, 2:end);
                w_proj_nTin = w_proj(~isTin_n);
                w_proj_nTin = [0 w_proj_nTin];
                
                % Weights without Tin & Tout at 450ms
                w_proj = w_dec{e, sess}(i_proj, 2:end);
                w_proj_nTiTo = w_proj(~isTiTo_n);
                w_proj_nTiTo = [0 w_proj_nTiTo];
                
                % Weights of when-decoder
                w_proj_whenC = [0 allW{7}{1, sess}(goodN, 1)'];
                w_proj_whenC_dim = [0 1];
                
                % Weights of when-decoder with just Tin neurons
                w_proj_whenC_Tin = [0 allW{70}{1, sess}(goodN, 1)'];
                
                % Weights of when-decoder without Tin neurons
                w_proj_whenC_nTin = [0 allW{71}{1, sess}(goodN, 1)'];
                
                for i = 1 : length(winStart)
                    
                    % Weights without Tin at same time point
                    w_proj = w_dec{e, sess}(i, 2:end);
                    w_proj_t = w_proj(~isTin_n);
                    w_proj_t = [0 w_proj_t];
                    
                    dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
                    
                    % Without Tin
                    xDD = squeeze(nanmean(spkZ.SL(useN, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    %xD = xDD(:, train&lDD)' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    % Weights without Tin at 450ms
                    choicePredTest  = glmval(w_proj_nTin',  -xTest, 'logit') > 0.5;
                    pc_nTin_450_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_nTin_450_test_se(i) = sqrt( (pc_nTin_450_test(i)*(1-pc_nTin_450_test(i))) / sum(~train&lDD)) ;
                    
                    % Weights without Tin at same time point
                    choicePredTest  = glmval(w_proj_t',  -xTest, 'logit') > 0.5;
                    pc_nTin_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_nTin_test_se(i) = sqrt( (pc_nTin_test(i)*(1-pc_nTin_test(i))) / sum(~train&lDD)) ;
                    
                    % ---------------------------------------------------
                    % Weights without Tin & Tout at same time point
                    w_proj = w_dec{e, sess}(i, 2:end);
                    w_proj_t = w_proj(~isTiTo_n);
                    w_proj_t = [0 w_proj_t];
                    
                    % Without Tin % Tout
                    xDD = squeeze(nanmean(spkZ.SL(useN_tito, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    % Weights without Tin at 450ms
                    choicePredTest  = glmval(w_proj_nTiTo',  -xTest, 'logit') > 0.5;
                    pc_nTiTo_450_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_nTiTo_450_test_se(i) = sqrt( (pc_nTiTo_450_test(i)*(1-pc_nTiTo_450_test(i))) / sum(~train&lDD)) ;
                    
                    % Weights without Tin at same time point
                    choicePredTest  = glmval(w_proj_t',  -xTest, 'logit') > 0.5;
                    pc_nTiTo_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_nTiTo_test_se(i) = sqrt( (pc_nTiTo_test(i)*(1-pc_nTiTo_test(i))) / sum(~train&lDD)) ;
                    
                    
                    % ---------------------------------------------------
                    % With Tin
                    xDD = squeeze(nanmean(spkZ.SL(goodN, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    %xD = xDD(:, train&lDD)' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    % Weights of when-decoder
                    choicePredTest  = glmval(w_proj_whenC',  -xTest, 'logit') > 0.5; % When decoder weights do not get applied to Z but 0 data!!!
                    pc_whenC_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_whenC_test_se(i) = sqrt( (pc_whenC_test(i)*(1-pc_whenC_test(i))) / sum(~train&lDD)) ;
                    % Weights of when-decoder using dimension traces
                    xDD_dim = nanmean(allDim_sl{7}{1, sess}(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); 
                    tTrain_dim = xDD_dim(train&lDD)' ; tTrain_dim = [zeros(size(tTrain_dim)); tTrain_dim];
                    [wDD,sDD]     =  lassoglm(-tTrain_dim', c', 'binomial','link','logit','Lambda',0);
                    tTest_dim = xDD_dim(~train&lDD)' ;
                    choicePredTest  = glmval(wDD,  -tTest_dim, 'logit') > 0.5; % When decoder weights do not get applied to Z but 0 data!!!
                    pc_whenC_dim_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_whenC_dim_test_se(i) = sqrt( (pc_whenC_dim_test(i)*(1-pc_whenC_dim_test(i))) / sum(~train&lDD)) ;
                    % When-decoder sign test
                    xDD_dim = nanmean(allDim_sl{7}{1, sess}(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); %lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    cho_sign = sign(xDD_dim); 
                    choicePredSign = cho_sign(lDD) == -1;
                    pc_whenC_sign_test(i)  =  sum(choicePredSign'==yChoice(lDD)) ./ sum(lDD) ;
                    
                    
                    % Weights of when-decoder with just Tin neurons
                    choicePredTest  = glmval(w_proj_whenC_Tin',  -xTest, 'logit') > 0.5;
                    pc_whenC_Tin_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_whenC_Tin_test_se(i) = sqrt( (pc_whenC_Tin_test(i)*(1-pc_whenC_Tin_test(i))) / sum(~train&lDD)) ;
                    % Weights of when-decoder using dimension traces
                    xDD_dim = nanmean(allDim_sl{70}{1, sess}(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); 
                    tTrain_dim = xDD_dim(train&lDD)' ; tTrain_dim = [zeros(size(tTrain_dim)); tTrain_dim];
                    [wDD,sDD]     =  lassoglm(-tTrain_dim', c', 'binomial','link','logit','Lambda',0);
                    tTest_dim = xDD_dim(~train&lDD)' ;
                    choicePredTest  = glmval(wDD,  -tTest_dim, 'logit') > 0.5; % When decoder weights do not get applied to Z but 0 data!!!
                    pc_whenC_Tin_dim_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_whenC_Tin_dim_test_se(i) = sqrt( (pc_whenC_Tin_dim_test(i)*(1-pc_whenC_Tin_dim_test(i))) / sum(~train&lDD)) ;
                    
                    % Weights of when-decoder without Tin neurons
                    choicePredTest  = glmval(w_proj_whenC_nTin',  -xTest, 'logit') > 0.5;
                    pc_whenC_nTin_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_whenC_nTin_test_se(i) = sqrt( (pc_whenC_nTin_test(i)*(1-pc_whenC_nTin_test(i))) / sum(~train&lDD)) ;
                    % Weights of when-decoder using dimension traces
                    xDD_dim = nanmean(allDim_sl{71}{1, sess}(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); 
                    tTrain_dim = xDD_dim(train&lDD)' ; tTrain_dim = [zeros(size(tTrain_dim)); tTrain_dim];
                    [wDD,sDD]     =  lassoglm(-tTrain_dim', c', 'binomial','link','logit','Lambda',0);
                    tTest_dim = xDD_dim(~train&lDD)' ;
                    choicePredTest  = glmval(wDD,  -tTest_dim, 'logit') > 0.5; % When decoder weights do not get applied to Z but 0 data!!!
                    pc_whenC_nTin_dim_test(i)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pc_whenC_nTin_dim_test_se(i) = sqrt( (pc_whenC_nTin_dim_test(i)*(1-pc_whenC_nTin_dim_test(i))) / sum(~train&lDD)) ;
                end
                pc450_Tin0_dec{e, sess} = pc_nTin_450_test;
                se450_Tin0_dec{e, sess} = pc_nTin_450_test_se;
                
                pc_Tin0_dec{e, sess} = pc_nTin_test;
                se_Tin0_dec{e, sess} = pc_nTin_test_se;
                
                pc450_TiTo0_dec{e, sess} = pc_nTiTo_450_test;
                se450_TiTo0_dec{e, sess} = pc_nTiTo_450_test_se;
                
                pc_TiTo0_dec{e, sess} = pc_nTiTo_test;
                se_TiTo0_dec{e, sess} = pc_nTiTo_test_se;
                
                pc_whenC_dec{e, sess} = pc_whenC_test;
                se_whenC_dec{e, sess} = pc_whenC_test_se;
                pc_whenC_sign_dec{e, sess} = pc_whenC_sign_test;
                
                pc_whenC_Tin_dec{e, sess} = pc_whenC_Tin_test;
                se_whenC_Tin_dec{e, sess} = pc_whenC_Tin_test_se;
                
                pc_whenC_nTin_dec{e, sess} = pc_whenC_nTin_test;
                se_whenC_nTin_dec{e, sess} = pc_whenC_nTin_test_se;
                
                pc_whenC_dim_dec{e, sess} = pc_whenC_dim_test;
                se_whenC_dim_dec{e, sess} = pc_whenC_dim_test_se;
                
                pc_whenC_Tin_dim_dec{e, sess} = pc_whenC_Tin_dim_test;
                se_whenC_Tin_dim_dec{e, sess} = pc_whenC_Tin_dim_test_se;
                
                pc_whenC_nTin_dim_dec{e, sess} = pc_whenC_nTin_dim_test;
                se_whenC_nTin_dim_dec{e, sess} = pc_whenC_nTin_dim_test_se;
            end
        end
    end
end
% The dim versions are needed for the when decoder
pc_whenC_dec = pc_whenC_dim_dec; 
se_whenC_dec = se_whenC_dim_dec; 
pc_whenC_Tin_dec = pc_whenC_Tin_dim_dec; 
se_whenC_Tin_dec = se_whenC_Tin_dim_dec; 
pc_whenC_nTin_dec = pc_whenC_nTin_dim_dec;
se_whenC_nTin_dec = se_whenC_nTin_dim_dec;


for sess = 1 : 8
    par = par_sess{sess};
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt>600; % l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = par.beh.cho_trg;
    for i = 1:length(winStart)
        dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
        % When-decoder sign test
        xDD_dim = nanmean(allDim_sl{7}{1, sess}(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); lDD = (~isnan(xDD_dim) == 1)';
        cho_sign = sign(xDD_dim);
        choicePredSign = cho_sign(lDD) == -1;
        pc_whenC_sign_test(i)  =  sum(choicePredSign==yChoice(lDD)) ./ sum(lDD) ;
    end
    pc_whenC_sign_dec{e, sess} = pc_whenC_sign_test;
end

save(fullfile(par_sess{1}.saveLoc, 'decoder_excluded_450_when'), 'w_dec', 'pc_dec', 'se_dec', 'excl', 'exclName', 'winCenter',...
    'pc450_dec', 'se450_dec', 'pc450_Tin0_dec', 'se450_Tin0_dec', 'pc_Tin0_dec', 'se_Tin0_dec',...
    'pc450_TiTo0_dec', 'se450_TiTo0_dec', 'pc_TiTo0_dec', 'se_TiTo0_dec',...
    'pc_whenC_dec', 'se_whenC_dec', 'pc_whenC_Tin_dec', 'se_whenC_Tin_dec', 'pc_whenC_nTin_dec', 'se_whenC_nTin_dec', 'pc_whenC_sign_dec')
% Should be same as 'E:\Matalb analyses\decoder_excluded_450_when_230522'



% figure; hold on; e = 1;
% for sess = 1 : 8
%     subplot(2, 4, sess); hold on
%     plot(winStart, pc_dec{e, sess}, 'k' , 'LineWidth', 2)
%     
%     plot(winStart, pc_whenC_dec{e, sess}, 'LineWidth', 2)
%     plot(winStart, pc_whenC_dim_dec{e, sess}, 'LineWidth', 2)
%     
%     plot(winStart, pc_whenC_Tin_dec{e, sess}, '--', 'LineWidth', 2)
%     plot(winStart, pc_whenC_Tin_dim_dec{e, sess}, '--', 'LineWidth', 2)
%     
%     plot(winStart, pc_whenC_nTin_dec{e, sess}, ':',  'LineWidth', 2)
%     plot(winStart, pc_whenC_nTin_dim_dec{e, sess}, ':', 'LineWidth', 2)
%     xlim([0 0.5])
%     ylim([0.4 0.85])
% end
% legend('all moving', 'whenC', 'whenC dim', 'whenC Tin', 'whenC Tin dim', 'whenC noTin', 'whenC noTin dim')


%% Grab weights at 450ms and make them into dimensions - response aligned

% load('E:\Matalb analyses\popMatZ_sessRL_short', 'rl_sessMatZ_short', 'trl_short')

i_proj = find(round(winCenter, 3) == 0.451); % i_proj = find(round(winCenter, 3) == 0.476); 
% i_proj = find(round(winCenter, 3) == 0.5510); 
% trl_short
trl_dec = -0.595 : 0.005 : 0; wBin = 0.005;

for sess = 1 : 8
    clear useDat
    for tp = 1 : length(trl_dec)
        tpts = find(trl_short >= trl_dec(tp) - wBin/2 & trl_short < trl_dec(tp) + wBin/2);
        useDat(:, :, tp) = nanmean(rl_sessMatZ_short{sess}(:, :, tpts), 3);
    end
    
    for e = 1 : size(pc_dec, 1)
        
        if size(w_dec{e, sess}, 1) >0
            w = w_dec{e, sess}(i_proj, 2:end); % used to be: w = w_dec{e, sess}(i_proj, :);
            %goodN = ones(size(rl_sessMatZ_short{sess}, 1), 1); goodN(excl{e}) = 0; goodN(sessZ{sess}.s==0) = 0;
            goodN = useN_dec{e, sess};
            
            C = squeeze(mmx('mult', w, useDat(goodN, :, :)));
            
            rl_exDec{e, sess} = C;
        end
    end
end


%% Response-aligned decoding performance!
%% 1. weights from 450ms applied to all time points
%% 2. weights from 450ms with Tin neurons set to zero applied to all time points
%% 3. weights with Tin neurons set to zero applied to same point
%% 4. weights of when decoder applied to all time points
%% 5. weights of when decoder with just Tin neurons applied to all time points
%% 6. weights of when decoder without Tin neurons applied to all time points


clear w_decR
for sess = 1 : 8
    par = par_sess{sess};
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = par.dataMat(l, par.idx.cho_trg);
    
    for e = 1 : length(excl)
        
        disp(num2str(sess))
        
        clear pcDD_test pDDtest_se useN
        
        useN = ones(size(rl_sessMatZ_short{sess}, 1), 1); useN(excl_sess{sess}{e}) = 0; useN(sessZ{sess}.s==0) = 0;
        useN_decR{e, sess} = find(useN==1);
        if sum(useN) > 0
            w_proj = w_dec{e, sess}(i_proj, :);
            
            for i = 1 : length(trl_dec)
                
                dW = find(trl_short >= trl_dec(i)-winSize/2 & trl_short <= trl_dec(i)+winSize/2); % Decoding window relative to dots onset
                
                xDD = squeeze(nanmean(rl_sessMatZ_short{sess}(useN==1, l, dW), 3));
                
                xDD = xDD - nanmean(xDD,2) ;
                lDD = (nanmean(~isnan(xDD), 1) == 1)';
                
                train = ~logical(mod(1:nL,2))' ;
                xD = xDD(:, train&lDD)' ;
                xTest = xDD(:, ~train&lDD)' ;
                c = yChoice(train&lDD);
                
                [wDD,sDD]     =  lassoglm(-xD, c, 'binomial','link','logit','Lambda',0.03);
                wDD = [0 ; wDD] ;
                
                choicePredTrain  = glmval(wDD,  -xD, 'logit') > 0.5;
                
                pcDD_train  =  sum(choicePredTrain==yChoice(train&lDD)) ./ sum(train&lDD) ;
                
                choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
                
                pcDD_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                pDDtest_se(i) = sqrt( (pcDD_test(i)*(1-pcDD_test(i))) / sum(~train&lDD)) ;
                
                w_decR{e, sess}(i, :) = wDD;
                
                choicePredTest  = glmval(w_proj',  -xTest, 'logit') > 0.5;
                
                pcDD550_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                pDDtest550_se(i) = sqrt( (pcDD550_test(i)*(1-pcDD550_test(i))) / sum(~train&lDD)) ;
                
            end
            pc_decR{e, sess} = pcDD_test;
            se_decR{e, sess} = pDDtest_se;
            
            pc450_decR{e, sess} = pcDD550_test;
            se450_decR{e, sess} = pDDtest550_se;
            
            if e == 1
                % Weights without Tin at same time point and weight without Tin from 450ms applied at each time point
                goodN = useN_dec{e, sess};
                isTin = zeros(1, size(sl_sessMatZ{sess}, 1)); isTin(par.unitIdxLIP_Tin) = 1; isTin_n = isTin(goodN);
                useN = goodN(~isTin_n);
                isTiTo = zeros(1, size(sl_sessMatZ{sess}, 1)); isTiTo([par.unitIdxLIP_Tin par.unitIdxLIP_Tout]) = 1; isTiTo_n = isTiTo(goodN);
                useN_tito = goodN(~isTiTo_n);
                
                % Weights without Tin at 450ms
                w_proj = w_dec{e, sess}(i_proj, 2:end);
                w_proj_nTin = w_proj(~isTin_n);
                w_proj_nTin = [0 w_proj_nTin];
                
                % Weights without Tin & Tout at 450ms
                w_proj = w_dec{e, sess}(i_proj, 2:end);
                w_proj_nTiTo = w_proj(~isTiTo_n);
                w_proj_nTiTo = [0 w_proj_nTiTo];
                
                %                 % Weights of when-decoder
                %                 w_proj_whenC = [0 allW{7}{1, sess}(goodN, 1)'];
                %
                %                 % Weights of when-decoder with just Tin neurons
                %                 w_proj_whenC_Tin = [0 allW{70}{1, sess}(goodN, 1)'];
                %
                %                 % Weights of when-decoder without Tin neurons
                %                 w_proj_whenC_nTin = [0 allW{71}{1, sess}(goodN, 1)'];
                
                
                for i = 1 : length(trl_dec)
                    
                    % Weights without Tin at same time point
                    w_proj = w_dec{e, sess}(i, 2:end);
                    w_proj_t = w_proj(~isTin_n);
                    w_proj_t = [0 w_proj_t];
                    
                    dW = find(trl_short >= trl_dec(i)-winSize/2 & trl_short <= trl_dec(i)+winSize/2);
                    dW_when = find(par.trl >= trl_dec(i)-winSize/2 & par.trl <= trl_dec(i)+winSize/2);
                    
                    % Without Tin
                    xDD = squeeze(nanmean(rl_sessMatZ_short{sess}(useN, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    %xD = xDD(:, train&lDD)' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    % Weights without Tin at 450ms
                    choicePredTest  = glmval(w_proj_nTin',  -xTest, 'logit') > 0.5;
                    pcR_nTin_450_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pcR_nTin_450_test_se(i) = sqrt( (pcR_nTin_450_test(i)*(1-pcR_nTin_450_test(i))) / sum(~train&lDD)) ;
                    
                    % Weights without Tin at same time point
                    choicePredTest  = glmval(w_proj_t',  -xTest, 'logit') > 0.5;
                    pcR_nTin_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pcR_nTin_test_se(i) = sqrt( (pcR_nTin_test(i)*(1-pcR_nTin_test(i))) / sum(~train&lDD)) ;
                    
                    % ---------------------------------------------------
                    % Weights without Tin & Tout at same time point
                    w_proj = w_dec{e, sess}(i, 2:end);
                    w_proj_t = w_proj(~isTiTo_n);
                    w_proj_t = [0 w_proj_t];
                    
                    % Without Tin % Tout
                    xDD = squeeze(nanmean(rl_sessMatZ_short{sess}(useN_tito, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
                    % Weights without Tin at 450ms
                    choicePredTest  = glmval(w_proj_nTiTo',  -xTest, 'logit') > 0.5;
                    pcR_nTiTo_450_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pcR_nTiTo_450_test_se(i) = sqrt( (pcR_nTiTo_450_test(i)*(1-pcR_nTiTo_450_test(i))) / sum(~train&lDD)) ;
                    
                    % Weights without Tin at same time point
                    choicePredTest  = glmval(w_proj_t',  -xTest, 'logit') > 0.5;
                    pcR_nTiTo_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                    pcR_nTiTo_test_se(i) = sqrt( (pcR_nTiTo_test(i)*(1-pcR_nTiTo_test(i))) / sum(~train&lDD)) ;
                    
                    
                    % ---------------------------------------------------
                    % With Tin
                    xDD = squeeze(nanmean(rl_sessMatZ_short{sess}(goodN, l, dW), 3));
                    xDD = xDD - nanmean(xDD,2) ;
                    lDD = (nanmean(~isnan(xDD), 1) == 1)';
                    
                    train = ~logical(mod(1:nL,2))' ;
                    %xD = xDD(:, train&lDD)' ;
                    xTest = xDD(:, ~train&lDD)' ;
                    c = yChoice(train&lDD);
                    
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
                pc450_Tin0_decR{e, sess} = pcR_nTin_450_test;
                se450_Tin0_decR{e, sess} = pcR_nTin_450_test_se;
                
                pc_Tin0_decR{e, sess} = pcR_nTin_test;
                se_Tin0_decR{e, sess} = pcR_nTin_test_se;
                
                pc450_TiTo0_decR{e, sess} = pcR_nTiTo_450_test;
                se450_TiTo0_decR{e, sess} = pcR_nTiTo_450_test_se;
                
                pc_TiTo0_decR{e, sess} = pcR_nTiTo_test;
                se_TiTo0_decR{e, sess} = pcR_nTiTo_test_se;
                
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

for sess = 1 : 8
    par = par_sess{sess};
    
    rt = par.beh.rt);
    coh = par.beh.dot_coh);
    l = rt>600; % l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = par.dataMat(l, par.idx.cho_trg);
    for i = 1:length(trl_dec)
        
        dW = find(par.trl >= trl_dec(i)-winSize/2 & par.trl <= trl_dec(i)+winSize/2);
        % When-decoder sign test
        xDD_dim = nanmean(allDim_rl{7}{1, sess}(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); lDD = (~isnan(xDD_dim) == 1)';
        cho_sign = sign(xDD_dim);
        choicePredSign = cho_sign(lDD) == -1;
        pcR_whenC_sign_test(i)  =  sum(choicePredSign==yChoice(lDD)) ./ sum(lDD) ;
    end
    pc_whenC_sign_decR{e, sess} = pcR_whenC_sign_test;
end


save('E:\Matalb analyses\decoder_excluded_RL_230522', 'w_dec', 'pc_dec', 'se_dec',...
    'pc_decR', 'se_decR', 'pc450_decR', 'se450_decR', 'excl', 'exclName',...
    'pc450_Tin0_decR', 'se450_Tin0_decR', 'pc_Tin0_decR', 'se_Tin0_decR',...
    'pc450_TiTo0_decR', 'se450_TiTo0_decR', 'pc_TiTo0_decR', 'se_TiTo0_decR',...
    'pc_whenC_decR', 'se_whenC_decR', 'pc_whenC_Tin_decR', 'se_whenC_Tin_decR', 'pc_whenC_nTin_decR', 'se_whenC_nTin_decR', 'pc_whenC_sign_decR')
% OLD
% save('E:\Matalb analyses\decoder_excluded_RL_230419', 'w_dec', 'pc_dec', 'se_dec',...
%     'pc_decR', 'se_decR', 'pc450_decR', 'se450_decR', 'excl', 'exclName',...
%     'pc450_Tin0_decR', 'se450_Tin0_decR', 'pc_Tin0_decR', 'se_Tin0_decR',...
%     'pc450_TiTo0_decR', 'se450_TiTo0_decR', 'pc_TiTo0_decR', 'se_TiTo0_decR',...
%     'pc_whenC_decR', 'se_whenC_decR', 'pc_whenC_Tin_decR', 'se_whenC_Tin_decR', 'pc_whenC_nTin_decR', 'se_whenC_nTin_decR')
% save('E:\Matalb analyses\decoder_excluded_RL_230508', 'w_dec', 'pc_dec', 'se_dec',...
%     'pc_decR', 'se_decR', 'pc450_decR', 'se450_decR', 'excl', 'exclName',...
%     'pc450_Tin0_decR', 'se450_Tin0_decR', 'pc_Tin0_decR', 'se_Tin0_decR',...
%     'pc450_TiTo0_decR', 'se450_TiTo0_decR', 'pc_TiTo0_decR', 'se_TiTo0_decR',...
%     'pc_whenC_decR', 'se_whenC_decR', 'pc_whenC_Tin_decR', 'se_whenC_Tin_decR', 'pc_whenC_nTin_decR', 'se_whenC_nTin_decR')


%% ==========================================================================
%% PLOTTING
%% Plot stimulus- and response-aligned traces for the same time point - Fig 1d, a and Supp Fig no-Tin
whenA = load('E:\Matalb analyses\choice_pred_from_Swhen');

rows = 1; cols = 2; 
xLimsSL = [0 0.6]; xLimsSL = [0 0.5]; 
xLimsRL = [-0.3 0]; 
ylims = [0.45 1];
pltC = {[1 9 14 12 13 2 16 5 15],... % [1 2 4 5 9 11 12 13],...
    [1 2 3 4 5],...
    };
sWin = 3;

% plotCollies{1} = {[0 0 0], [55 158 176]/255, [], [76 0 153]/255,...
%     [255 102 102]/255, [0 153 76]/255, [20 20 255]/255,...
%     [204 102 0]/255, [55 158 176]/255, [255 102 102]/255, [102 0 204]/255,...
%     [255 102 102]/255, [153 0 76]/255, [0 153 76]/255, [20 20 255]/255,...
%     [204 102 0]/255,...
%     };% [153 0 76]/255, [255 200 0]/255,
plotCollies{1} = {[0 0 0], [55 158 176]/255, [], [76 0 153]/255,...
    [255 145 105]/255, [0 153 76]/255, [20 20 255]/255,...
    [204 102 0]/255, [55 158 176]/255, [255 102 102]/255, [102 0 204]/255,...
    [255 145 105]/255,...
    [200 50 53]/255, [76 0 153]/255, [200 50 53]/255, [76 0 153]/255,... % justTinMin, justMinCI, noTinMin, noMinCI
    [0 153 76]/255,...
    [20 20 255]/255,...
    [204 102 0]/255,...
    };% [153 0 76]/255, [255 200 0]/255, [153 0 76]/255, [255 102 102]/255

% [200 50 53]
% [255 175 143]
plotCollies{2} = plotCollies{1};
plotCollies{2}{1} = [176 28 98]/255;
d = [1 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2]; dashes = {'-', ':'};
lws = [3 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1];
lws2 = [3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

clear pltVar
for sess = 1 : 8
    v = 0;
    v=v+1;pltVar{v}.m(:, sess) = pc_dec{1, sess};
    pltVar{v}.s(:, sess) = se_dec{1, sess};
    pltVar{v}.mr(:, sess) = pc_decR{1, sess};
    pltVar{v}.sr(:, sess) = se_decR{1, sess};
    pltNames{v} = 'all moving window';
    pltCols{v} = [243 158 17]/255;
    pltDash{v} = '-';
    
    % Figure 4a
    v=v+1;pltVar{v}.m(:, sess) = pc450_dec{1, sess};
    pltVar{v}.s(:, sess) = se450_dec{1, sess};
    pltVar{v}.mr(:, sess) = pc450_decR{1, sess};
    pltVar{v}.sr(:, sess) = se450_decR{1, sess};
    pltNames{v} = 'all fixed window';
    pltCols{v} = [176 28 98]/255;
    pltDash{v} = '-';
    
    %     v=v+1;pltVar{v}.m(:, sess) = pc_TiTo0_dec{1, sess};
    %     pltVar{v}.s(:, sess) = se_TiTo0_dec{1, sess};
    %     pltVar{v}.mr(:, sess) = pc_TiTo0_decR{1, sess};
    %     pltVar{v}.sr(:, sess) = se_TiTo0_decR{1, sess};
    %     pltNames{v} = 'T_i_n & T_o_u_t zero';
    %     pltCols{v} = [255 145 105]/255;
    %     pltDash{v} = '-';
    
    %     v = v + 1;
    %     dtA = whenA.acc_pred.motion_aligned.t(2) - whenA.acc_pred.motion_aligned.t(1);
    %     dtW = winCenter(2) - winCenter(1);
    %     pltDat = nan(size(winCenter));
    %     useTPT = [find(whenA.acc_pred.motion_aligned.t == winCenter(1)) : round(dtW/dtA) : whenA.acc_pred.motion_aligned.t(end)/dtA];
    %     pltDat(1 : length(useTPT)) = whenA.acc_pred.motion_aligned.pred(useTPT);
    %     pltVar{v}.m(:, sess) = pltDat;
    %     pltVar{v}.s(:, sess) = nan(size(pltDat));
    %
    %     pltDat = nan(size(trl_dec));
    %     useTPT = [round(dtW/dtA) : round(dtW/dtA) : find(whenA.acc_pred.resp_aligned.t == trl_dec(end))];
    %     pltDat(end-length(useTPT)+1:end) = whenA.acc_pred.resp_aligned.pred(useTPT);
    %     pltVar{v}.mr(:, sess) = pltDat;
    %     pltVar{v}.sr(:, sess) = nan(size(pltDat));
    %     pltNames{v} = 'when decoder (sign)';
    %     pltCols{v} = [41 174 137]/255;
    %     pltDash{v} = '-';
    
    v=v+1;pltVar{v}.m(:, sess) = pc_whenC_sign_dec{1, sess};
    pltVar{v}.s(:, sess) = nan(size(se_whenC_dec{1, sess}));
    pltVar{v}.mr(:, sess) = pc_whenC_sign_decR{1, sess};
    pltVar{v}.sr(:, sess) = nan(size(se_whenC_dec{1, sess}));
    pltNames{v} = 'when decoder (sign)';
    pltCols{v} = [41 174 137]/255;
    pltDash{v} = '-';
    
    
    v=v+1;pltVar{v}.m(:, sess) = pc_whenC_dec{1, sess};
    pltVar{v}.s(:, sess) = se_whenC_dec{1, sess};
    pltVar{v}.mr(:, sess) = pc_whenC_decR{1, sess};
    pltVar{v}.sr(:, sess) = se_whenC_decR{1, sess};
    pltNames{v} = 'when decoder (regression)';
    pltCols{v} = [41 174 137]/255;
    pltDash{v} = '--';
    
    % Figure Supp no-Tin
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{5, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{5, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{5, sess};
%     pltVar{v}.sr(:, sess) = se_decR{5, sess};
%     pltNames{v} = 'no T_i_n & T_o_u_t refit';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{5, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{5, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{5, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{5, sess};
%     pltNames{v} = 'no T_i_n & T_o_u_t refit fixed';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '--';
%     
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{12, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{12, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{12, sess};
%     pltVar{v}.sr(:, sess) = se_decR{12, sess};
%     pltNames{v} = 'only T_i_n & T_o_u_t refit';
%     pltCols{v} = [176 28 98]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{12, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{12, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{12, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{12, sess};
%     pltNames{v} = 'only T_i_n & T_o_u_t refit fixed';
%     pltCols{v} = [176 28 98]/255;
%     pltDash{v} = '--';
%     
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{14, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{14, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{14, sess};
%     pltVar{v}.sr(:, sess) = se_decR{14, sess};
%     pltNames{v} = 'no T_i_n & M_i_n';
%     pltCols{v} = [171 145 2]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{14, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{14, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{14, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{14, sess};
%     pltNames{v} = 'no T_i_n & M_i_n fixed';
%     pltCols{v} = [171 145 2]/255;
%     pltDash{v} = '--';
    
    
    
%     v=v+1;pltVar{v}.m(:, sess) = pc_Tin0_dec{1, sess};
%     pltVar{v}.s(:, sess) = se_Tin0_dec{1, sess};
%     pltVar{v}.mr(:, sess) = pc_Tin0_decR{1, sess};
%     pltVar{v}.sr(:, sess) = se_Tin0_decR{1, sess};
%     pltNames{v} = 'T_i_n zero';
%     pltCols{v} = [255 145 105]/255;
%     pltDash{v} = '-';
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{2, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{2, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{2, sess};
%     pltVar{v}.sr(:, sess) = se_decR{2, sess};
%     pltNames{v} = 'no T_i_n refit';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{2, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{2, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{2, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{2, sess};
%     pltNames{v} = 'no T_i_n refit fixed';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '--';
    
%     v=v+1;pltVar{v}.m(:, sess) = pc_whenC_dec{1, sess};
%     pltVar{v}.s(:, sess) = se_whenC_dec{1, sess};
%     pltVar{v}.mr(:, sess) = pc_whenC_decR{1, sess};
%     pltVar{v}.sr(:, sess) = se_whenC_decR{1, sess};
%     pltNames{v} = 'when decoder';
%     pltCols{v} = [237 190 4]/255;
%     pltDash{v} = '-';
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_whenC_Tin_dec{1, sess};
%     pltVar{v}.s(:, sess) = se_whenC_Tin_dec{1, sess};
%     pltVar{v}.mr(:, sess) = pc_whenC_Tin_decR{1, sess};
%     pltVar{v}.sr(:, sess) = se_whenC_Tin_decR{1, sess};
%     pltNames{v} = 'when T_i_n only';
%     pltCols{v} = [41 174 137]/255;
%     pltDash{v} = '-';
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_whenC_nTin_dec{1, sess};
%     pltVar{v}.s(:, sess) = se_whenC_nTin_dec{1, sess};
%     pltVar{v}.mr(:, sess) = pc_whenC_nTin_decR{1, sess};
%     pltVar{v}.sr(:, sess) = se_whenC_nTin_decR{1, sess};
%     pltNames{v} = 'when no T_i_n';
%     pltCols{v} = [76 0 153]/255;
%     pltDash{v} = '-';
    
end

% d = ;

figure('Position', [600 400 700 300]); hold on;
% -------------------------------
% Stimulus-locked - same time point
sfh1 = subplot(rows, cols, 1); hold on
legs = []; l = 0;
for pl = 1 : length(pltVar)
    pltMat = pltVar{pl}.m;
    plot(winCenter, smooth(nanmean(pltMat, 2), sWin),...
        'Color', pltCols{pl}, 'LineWidth', 1, 'LineStyle', pltDash{pl}) % lws2(e)
    l = l + 1; legs{l} = pltNames{pl};
end
title(['mean across sessions'])
xlabel('time to response')
ylabel('xval decoder perf')
xlim(xLimsSL)
ylim(ylims)
legend(legs, 'Location', 'NorthWest')
grid on
sfh1.Position = [0.1300 0.15 0.35 0.775];
slWid = sfh1.Position(3);
% ====================================
% Response-locked - same time point
sfh2 = subplot(rows, cols, 2); hold on
legs = []; l = 0;
legs = []; l = 0;
for pl = 1 : length(pltVar)
    pltMat = pltVar{pl}.mr;
    plot(trl_dec, smooth(nanmean(pltMat, 2), sWin),...
        'Color', pltCols{pl}, 'LineWidth', 1, 'LineStyle', pltDash{pl}) % lws2(e)
    l = l + 1; legs{l} = pltNames{pl};
end
title(['mean across sessions'])
xlabel('time to response')
% ylabel('xval decoder perf')
xlim(xLimsRL)
ylim(ylims)
% legend(legs, 'Location', 'NorthWest')
grid on
rlWid = slWid / abs(diff(xLimsSL)) * abs(diff(xLimsRL));
sfh2.Position = [0.5703 0.15 rlWid 0.775];



%% Plot stimulus and response aligned decoding accuracy for the when decoder with and without Tins

clear pltVar
for sess = 1 : 8
    v = 0;
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{1, sess};
%     pltVar{v}.s(:, sess) = se_dec{1, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{1, sess};
%     pltVar{v}.sr(:, sess) = se_decR{1, sess};
%     pltNames{v} = 'all moving window';
%     pltCols{v} = 0.4*[1 1 1];
%     pltDash{v} = '-';
    
    v=v+1;pltVar{v}.m(:, sess) = pc450_dec{1, sess};
    pltVar{v}.s(:, sess) = se450_dec{1, sess};
    pltVar{v}.mr(:, sess) = pc450_decR{1, sess};
    pltVar{v}.sr(:, sess) = se450_decR{1, sess};
    pltNames{v} = 'all fixed window';
    pltCols{v} = 0.2*[1 1 1];
    pltDash{v} = '-';
     
%     v=v+1;pltVar{v}.m(:, sess) = pc_TiTo0_dec{1, sess};
%     pltVar{v}.s(:, sess) = se_TiTo0_dec{1, sess};
%     pltVar{v}.mr(:, sess) = pc_TiTo0_decR{1, sess};
%     pltVar{v}.sr(:, sess) = se_TiTo0_decR{1, sess};
%     pltNames{v} = 'T_i_n & T_o_u_t zero';
%     pltCols{v} = [255 145 105]/255;
%     pltDash{v} = '-';
    
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{5, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{5, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{5, sess};
%     pltVar{v}.sr(:, sess) = se_decR{5, sess};
%     pltNames{v} = 'no T_i_n & T_o_u_t refit';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{5, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{5, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{5, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{5, sess};
%     pltNames{v} = 'no T_i_n & T_o_u_t refit fixed';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '--';
%     
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{12, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{12, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{12, sess};
%     pltVar{v}.sr(:, sess) = se_decR{12, sess};
%     pltNames{v} = 'only T_i_n & T_o_u_t refit';
%     pltCols{v} = [176 28 98]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{12, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{12, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{12, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{12, sess};
%     pltNames{v} = 'only T_i_n & T_o_u_t refit fixed';
%     pltCols{v} = [176 28 98]/255;
%     pltDash{v} = '--';
%     
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{14, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{14, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{14, sess};
%     pltVar{v}.sr(:, sess) = se_decR{14, sess};
%     pltNames{v} = 'no T_i_n & M_i_n';
%     pltCols{v} = [171 145 2]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{14, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{14, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{14, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{14, sess};
%     pltNames{v} = 'no T_i_n & M_i_n fixed';
%     pltCols{v} = [171 145 2]/255;
%     pltDash{v} = '--';
    
    
    
%     v=v+1;pltVar{v}.m(:, sess) = pc_Tin0_dec{1, sess};
%     pltVar{v}.s(:, sess) = se_Tin0_dec{1, sess};
%     pltVar{v}.mr(:, sess) = pc_Tin0_decR{1, sess};
%     pltVar{v}.sr(:, sess) = se_Tin0_decR{1, sess};
%     pltNames{v} = 'T_i_n zero';
%     pltCols{v} = [255 145 105]/255;
%     pltDash{v} = '-';
%     
%     v=v+1;pltVar{v}.m(:, sess) = pc_dec{2, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se_dec{2, sess};
%     pltVar{v}.mr(:, sess) = pc_decR{2, sess};
%     pltVar{v}.sr(:, sess) = se_decR{2, sess};
%     pltNames{v} = 'no T_i_n refit';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '-';
%     
%     v=v+1; pltVar{v}.m(:, sess) = pc450_dec{2, sess}; % no Tin refit
%     pltVar{v}.s(:, sess) = se450_dec{2, sess};
%     pltVar{v}.mr(:, sess) = pc450_decR{2, sess};
%     pltVar{v}.sr(:, sess) = se450_decR{2, sess};
%     pltNames{v} = 'no T_i_n refit fixed';
%     pltCols{v} = [35 108 244]/255;
%     pltDash{v} = '--';
    
    v=v+1;pltVar{v}.m(:, sess) = pc_whenC_dec{1, sess};
    pltVar{v}.s(:, sess) = se_whenC_dec{1, sess};
    pltVar{v}.mr(:, sess) = pc_whenC_decR{1, sess};
    pltVar{v}.sr(:, sess) = se_whenC_decR{1, sess};
    pltNames{v} = 'when decoder';
    pltCols{v} = [237 190 4]/255;
    pltDash{v} = '-';
    
    v=v+1;pltVar{v}.m(:, sess) = pc_whenC_Tin_dec{1, sess};
    pltVar{v}.s(:, sess) = se_whenC_Tin_dec{1, sess};
    pltVar{v}.mr(:, sess) = pc_whenC_Tin_decR{1, sess};
    pltVar{v}.sr(:, sess) = se_whenC_Tin_decR{1, sess};
    pltNames{v} = 'when T_i_n only';
    pltCols{v} = [41 174 137]/255;
    pltDash{v} = '-';
    
    v=v+1;pltVar{v}.m(:, sess) = pc_whenC_nTin_dec{1, sess};
    pltVar{v}.s(:, sess) = se_whenC_nTin_dec{1, sess};
    pltVar{v}.mr(:, sess) = pc_whenC_nTin_decR{1, sess};
    pltVar{v}.sr(:, sess) = se_whenC_nTin_decR{1, sess};
    pltNames{v} = 'when no T_i_n';
    pltCols{v} = [76 0 153]/255;
    pltDash{v} = '-';
    
end

% d = ;

figure('Position', [600 400 700 300]); hold on;
% -------------------------------
% Stimulus-locked - same time point
sfh1 = subplot(rows, cols, 1); hold on
legs = []; l = 0;
for pl = 1 : length(pltVar)
    pltMat = pltVar{pl}.m;
    plot(winCenter, smooth(nanmean(pltMat, 2), sWin),...
        'Color', pltCols{pl}, 'LineWidth', 1, 'LineStyle', pltDash{pl}) % lws2(e)
    l = l + 1; legs{l} = pltNames{pl};
end
title(['mean across sessions'])
xlabel('time to response')
ylabel('xval decoder perf')
xlim(xLimsSL)
ylim(ylims)
legend(legs, 'Location', 'NorthWest')
grid on
sfh1.Position = [0.1300 0.15 0.35 0.775];
slWid = sfh1.Position(3);
% ====================================
% Response-locked - same time point
sfh2 = subplot(rows, cols, 2); hold on
legs = []; l = 0;
for pl = 1 : length(pltVar)
    pltMat = pltVar{pl}.mr;
    plot(trl_dec, smooth(nanmean(pltMat, 2), sWin),...
        'Color', pltCols{pl}, 'LineWidth', 1, 'LineStyle', pltDash{pl}) % lws2(e)
    l = l + 1; legs{l} = pltNames{pl};
end
title(['mean across sessions'])
xlabel('time to response')
% ylabel('xval decoder perf')
xlim(xLimsRL)
ylim(ylims)
% legend(legs, 'Location', 'NorthWest')
grid on
rlWid = slWid / abs(diff(xLimsSL)) * abs(diff(xLimsRL));
sfh2.Position = [0.5703 0.15 rlWid 0.775];




%% Plot stimulus- and response-aligned traces for the same time point
rows = 1; cols = 2; 
xLimsSL = [0 0.6]; xLimsSL = [0 0.5]; 
xLimsRL = [-0.3 0]; 
ylims = [0.45 1];
pltC = {[1 9 14 12 13 2 16 5 15],... % [1 2 4 5 9 11 12 13],...
    [1 2 3 4 5],...
    };
sWin = 3;

plotCollies{1} = {[0 0 0], [55 158 176]/255, [], [76 0 153]/255,...
    [255 145 105]/255, [0 153 76]/255, [20 20 255]/255,...
    [204 102 0]/255, [55 158 176]/255, [255 102 102]/255, [102 0 204]/255,...
    [255 145 105]/255,...
    [200 50 53]/255, [76 0 153]/255, [200 50 53]/255, [76 0 153]/255,... % justTinMin, justMinCI, noTinMin, noMinCI
    [0 153 76]/255,...
    [20 20 255]/255,...
    [204 102 0]/255,...
    };

plotCollies{2} = plotCollies{1};
plotCollies{2}{1} = [176 28 98]/255;
d = [1 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2]; dashes = {'-', ':'};
lws = [3 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1];
lws2 = [3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

use450 = 1;

figure('Position', [600 200 700 300]); hold on;
% -------------------------------
% Stimulus-locked - same time point
sfh1 = subplot(rows, cols, 1); hold on
legs = []; l = 0;
c = 1;
for dd = 1 % : 2
    for e = pltC{c}
        pltMat = nan(8, size(pc_dec{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_dec{e, sess}, 2)>0
                if use450 == 0
                    pltMat(sess, :) = pc_dec{e, sess};
                else
                    pltMat(sess, :) = pc450_dec{e, sess};
                end
            end
        end
        plot(winCenter, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCollies{1}{e}, 'LineWidth', 1, 'LineStyle', dashes{d(e)}) % lws2(e)
        l = l + 1; legs{l} = exclName{e};
    end
end
title(['mean across sessions'])
xlabel('time to response')
ylabel('xval decoder perf')
xlim(xLimsSL)
ylim(ylims)
legend(legs, 'Location', 'NorthWest')
grid on
sfh1.Position = [0.1300 0.15 0.35 0.775];
slWid = sfh1.Position(3);
% ====================================
% Response-locked - same time point
sfh2 = subplot(rows, cols, 2); hold on
legs = []; l = 0;
c = 1;
for dd = 1 % : 2
    for e = pltC{c}
        pltMat = nan(8, size(pc_decR{1, 1}, 2));
        for sess = 1 : 8
            
            if size(pc_decR{e, sess}, 2)>0
                if use450 == 0
                    pltMat(sess, :) = pc_decR{e, sess};
                else
                    pltMat(sess, :) = pc450_decR{e, sess};
                end
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCollies{1}{e}, 'LineWidth', 1, 'LineStyle', dashes{d(e)}) % lws2(e)
        l = l + 1; legs{l} = exclName{e};
    end
end
title(['mean across sessions'])
xlabel('time to response')
% ylabel('xval decoder perf')
xlim(xLimsRL)
ylim(ylims)
% legend(legs, 'Location', 'NorthWest')
grid on
rlWid = slWid / abs(diff(xLimsSL)) * abs(diff(xLimsRL));
sfh2.Position = [0.5703 0.15 rlWid 0.775];

%% Bar graph comparing performance at -100ms between moving and fixed windows
useTP = find(trl_dec==-0.1);

clear barDat mov fix
for e = 1 : size(pc450_decR, 1)
    for sess = 1 : 8
        if size(pc450_decR{e, sess}, 1) > 0
            a = pc_decR{e, sess}(1, useTP);
            b = pc450_decR{e, sess}(1, useTP);
            
            barDat(e, sess) = (a - b) / a;
            mov(e, sess) = a;
            fix(e, sess) = b;
        else
            barDat(e, sess) = nan;
            mov(e, sess) = nan;
            fix(e, sess)  = nan;
        end
    end
end

pltE = [1 5 12 14];

figure; hold on
i = 0; 
for p = pltE
    i = i + 1;
    m = nanmean(barDat(p, :), 2);
    s = nanstd(barDat(p, :)) / sqrt(sum(~isnan(barDat(p, :))));
    plot(i, m, 'o')
    line([i i], m + [-s s], 'Color', 'k')
end

pltCol = {131 * [1 1 1]/255,...
    [170 28 34]/255,...
    plotCollies{1}{12},...
    plotCollies{1}{14}};


figure; hold on
x = [1 : length(pltE)];
m = nanmean(barDat(pltE, :), 2);
s = nanstd(barDat(pltE, :), [], 2) ./ sqrt(8);
b = bar(x, m,'FaceColor','flat') ;
for p = 1 : length(pltE)
    b.CData(p,:) = pltCol{p};
end
er = errorbar(x,m,s,s);
er.Color = [0 0 0];
er.LineStyle = 'none';  
hold off
xlim([0.5 length(pltE)+0.5])

%% Compute the decoder performance of the when-decoder at the time of choice committment
% Computed differently
useDims = [7 9 70 71];
useT = find(par.trl >= -0.125 & par.trl <= -0.075);

d = 0;
for dim = useDims
    d = d + 1;
    for sess = 1 : 8
        par = par_sess{sess};
        
        rt = par.beh.rt;
        coh = par.beh.dot_coh;
        l = rt>600; % l = rt>600  & coh >= 0;
        nL = sum(l);
        yChoice = par.dataMat(l, par.idx.cho_trg);
        
        xDD = nanmean(allDim_rl{dim}{1, sess}(l, useT), 2)';
        
        xDD = xDD - nanmean(xDD,2) ;
        lDD = (nanmean(~isnan(xDD), 1) == 1)';
        
        train = ~logical(mod(1:nL,2))' ;
        xD = xDD(:, train&lDD)' ;
        xTest = xDD(:, ~train&lDD)' ;
        c = yChoice(train&lDD);
        
        [wDD,sDD]     =  lassoglm(-xD, c, 'binomial','link','logit','Lambda',0.03);
        wDD = [0 ; wDD] ;
        
        choicePredTrain  = glmval(wDD,  -xD, 'logit') > 0.5;
        pcDD_train  =  sum(choicePredTrain==yChoice(train&lDD)) ./ sum(train&lDD) ;
        
        choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
        pcDD_test  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
        
        decPerf(d, sess) = pcDD_test;
        
    end
end

for d = 1 : length(useDims)
    disp([num2str(dimName{useDims(d)}) ': ' num2str(round(nanmean(decPerf(d, :)), 2)) ' +- ' num2str(round(nanstd(decPerf(d, :))/sqrt(8), 2))])
end
% reg2whenC: 0.89 +- 0.03
% reg2whenC_Tin0: 0.79 +- 0.04
% reg2whenC_justTin: 0.86 +- 0.03
% reg2whenC_noTin: 0.81 +- 0.05

pCols = {0.5*[1 1 1], [41 174 137]/255, [237 190 4]/255};
plot_d = [1 3 4];
figure; hold on
bar([1 2 3],nanmean(decPerf(plot_d, :), 2));
for d = 1 : length(plot_d)
    hold on;
    m = nanmean(decPerf(plot_d(d), :), 2);
    sem = nanstd(decPerf(plot_d(d), :))/sqrt(8);
    bar(d,m,'FaceColor',pCols{d});
    line(d*[1 1], m+[-sem sem], 'Color', 'k')
end
xlim([0.5 length(plot_d)+0.5])
ylim([0.5 1])

%% Plot stimulus- and response-aligned traces for the same time point and weights derived at 450ms
rows = 1; cols = 2; 
xLimsSL = [0 0.6]; 
xLimsRL = [-0.3 0]; 
ylims = [0.45 1];
pltC = {[1 5 12 15],... % [1 2 5 9 11 12 13],...
    [1 2 3 4 5],...
    };
exclName = {'All', 'noTin', 'noTout', 'noMinC', 'noTiTo', 'noTiMi', 'noToMi', 'noTiToMi',...
    'justTin', 'justTout', 'justMin', 'justTiTo', 'noTinMin', 'noMinCI', 'justMinCI'}; 
plotCollies{1} = {[0 0 0], [55 158 176]/255, [], [76 0 153]/255,...
    [255 102 102]/255, [0 153 76]/255, [20 20 255]/255,...
    [204 102 0]/255, [55 158 176]/255, [255 102 102]/255, [102 0 204]/255,...
    [255 102 102]/255, [153 0 76]/255, [0 153 76]/255, [20 20 255]/255,...
    [204 102 0]/255,...
    };% [153 0 76]/255, [255 200 0]/255,
plotCollies{2} = plotCollies{1};
plotCollies{2}{1} = [176 28 98]/255;

figure('Position', [600 200 700 300]); hold on; 
% -------------------------------
% Stimulus-locked - same time point
sfh1 = subplot(rows, cols, 1); hold on
legs = []; l = 0;
c = 1;
for dd = 1 : 2
    for e = pltC{c}
        pltMat = nan(8, size(pc_dec{1, 1}, 2));
        for sess = 1 : 8
            if dd == 1
                if size(pc_dec{e, sess}, 2)>0
                    pltMat(sess, :) = pc_dec{e, sess};
                end
            elseif dd == 2
                if size(pc450_dec{e, sess}, 2)>0
                    pltMat(sess, :) = pc450_dec{e, sess};
                end
            end
        end
        plot(winCenter, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCollies{dd}{e}, 'LineWidth', 1, 'LineStyle', dashes{d(e)}) % lws2(e)
        l = l + 1; legs{l} = exclName{e};
    end
end
title(['mean across sessions'])
xlabel('time to response')
ylabel('xval decoder perf')
xlim(xLimsSL)
ylim(ylims)
legend('all moving', 'Tin moving', 'no Tin moving', 'no Min moving',...
    'all fixed', 'Tin fixed', 'no Tin fixed', 'no Min fixed',...
    'NorthWest')
%legend('moving window', 'fixed window', 'Location', 'NorthWest')
grid on
sfh1.Position = [0.1300 0.15 0.35 0.775];
slWid = sfh1.Position(3);
% ====================================
% Response-locked - same time point
sfh2 = subplot(rows, cols, 2); hold on
legs = []; l = 0;
c = 1;
for dd = 1 : 2
    for e = pltC{c}
        pltMat = nan(8, size(pc_decR{1, 1}, 2));
        for sess = 1 : 8
             if dd == 1
                if size(pc_decR{e, sess}, 2)>0
                    pltMat(sess, :) = pc_decR{e, sess};
                end
            elseif dd == 2
                if size(pc450_decR{e, sess}, 2)>0
                    pltMat(sess, :) = pc450_decR{e, sess};
                end
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCollies{dd}{e}, 'LineWidth', 1, 'LineStyle', dashes{d(e)}) % lws2(e)
        l = l + 1; legs{l} = exclName{e};
    end
end
title(['mean across sessions'])
xlabel('time to response')
% ylabel('xval decoder perf')
xlim(xLimsRL)
ylim(ylims)
% legend(legs, 'Location', 'NorthWest')
grid on
rlWid = slWid / abs(diff(xLimsSL)) * abs(diff(xLimsRL));
sfh2.Position = [0.5703 0.15 rlWid 0.775];



%% =======================================================================
% Run decoder on simulated data
% simDat = load('C:\Users\shadlenlab\Dropbox\Natalie_and_Mike\simulation4NSt'); 
simDat = load('C:\Users\shadlenlab\Dropbox\Natalie_and_Mike\sim4Natalie-29Apr'); 

delay = 0.2;
time = round(simDat.P.t+delay, 3);
minRT = 0.6;
tMin = find(time == minRT);
simDV = [];
simCho = [];
simRT = [];
for cc = 1 : 11
    newDV = squeeze(simDat.D.dv(:, cc, :));
    simDV = [simDV; newDV];
    simCho = [simCho simDat.D.ch(:, cc)'];
    simRT = [simRT simDat.D.decT(:, cc)'];
end
% useTrl = find(~isnan(simDV(:, tMin)));
useTrl = find(simRT+delay >= minRT);

for i = 1 : length(winStart)
    
    yChoice = simCho(useTrl)';
    dW = find(time >= winStart(i) & time <= winEnd(i)); % Decoding window relative to dots onset
    if length(dW) > 0
        xDD = nanmean(simDV(useTrl, dW), 2);
        
        xDD = xDD - nanmean(xDD) ;
        lDD = ~isnan(xDD);
        nL = sum(lDD);
        
%         train = ~logical(mod(1:nL,2))' ;
%         %xD = xDD(:, train&lDD)' ;
%         xTest = xDD(~train&lDD)' ;
%         c = yChoice(train&lDD);
        
        choicePredTest  = glmval(1,  xDD, 'logit') > 0.5;
        
        pcSim(i)  =  sum(choicePredTest(:, 2)==yChoice) ./ length(yChoice) ;
        pcSim_se(i) = sqrt( (pcSim(i)*(1-pcSim(i))) / length(yChoice)) ;
    else
        pcSim(i)  =  nan ;
        pcSim_se(i) = nan ;
    end
    
    % pcSim(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
    % pcSim_se(i) = sqrt( (pcSim(i)*(1-pcSim(i))) / sum(~train&lDD)) ;
end

tw = winStart >= delay;
figure; hold on
plot(winStart(tw), pcSim(tw))
xlim([0 0.5])
ylim([0.45 1])


figure; hold on
subplot(1, 3, 1); hold on
plot(winStart(tw), pcSim(tw))
xlim([0 0.3])
ylim([0.45 1])
subplot(1, 3, 2); hold on
plot(winStart(tw), pcSim(tw))
xlim([0 0.5])
ylim([0.45 1])
subplot(1, 3, 3); hold on
plot(winStart(tw), pcSim(tw))
xlim([0.2 0.5])
ylim([0.45 1])


save('E:\Matalb analyses\sl_decoder_simulated', 'pcSim', 'pcSim_se', 'simDV', 'simCho', 'useTrl', 'time')



%% =======================================================================
%% OTHER PLOTS
%% Plot example session and mean across sessions
sWin = 3; xlims = [0.1 0.5];
figure; hold on
subplot(1, 2, 1); hold on
for sess = 1
    e = 1; plot(winCenter, smooth(pc_dec{e, sess}, sWin), 'k', 'LineWidth', 2)
    for e = 2 : 8
        plot(winCenter, smooth(pc_dec{e, sess}, sWin))
    end
    for e = 9 : 13
        plot(winCenter, smooth(pc_dec{e, sess}, sWin), ':', 'LineWidth', 2)
    end
    % legend(exclName)
    title(par_sess{sess}.date)
    xlabel('time from motin onset')
    ylabel('xval decoder perf')
end
xlim(xlims)
%--
subplot(1, 2, 2); hold on
e = 1; pltMat = []; l = 0;
for sess = 1 : 8
    pltMat(sess, :) = pc_dec{e, sess};
end
plot(winCenter, nanmean(pltMat, 1), 'k', 'LineWidth', 2)
l = l + 1; legs{l} = exclName{e};
for e = 2 : 8
    for sess = 1 : 8
        pltMat(sess, :) = pc_dec{e, sess};
    end
    plot(winCenter, nanmean(pltMat, 1))
end
for e = 9 : 13
    pltMat = nan(8, size(pc_dec{1, 1}, 2));
    for sess = 1 : 8
        if length(pc_dec{e, sess}) > 0
            pltMat(sess, :) = pc_dec{e, sess};
        end
    end
    plot(winCenter, nanmean(pltMat, 1), ':', 'LineWidth', 2)
end
legend(exclName)
title('mean across sessions')
xlabel('time from motin onset')
ylabel('xval decoder perf')
xlim(xlims)
% ====================================
% Plot lines relevant to each class of neurons

rows = 3; cols = 2; 
sWin = 3; xlims = [0.1 0.5]; ylims = [0.45 0.75];
plotCols = {0.3*[1 1 1],...
    [20 20 255]/255, [20 20 255]/255, [20 20 255]/255,...
    [153 0 76]/255, [255 102 102]/255, [255 200 0]/255, [102 0 204]/255,...
    [0 153 76]/255, [0 153 76]/255, [0 153 76]/255,...
    [0 255 255]/255,...
    };
plotCols2 = {0.3*[1 1 1],...
    [153 0 76]/255, [0 153 76]/255, [20 20 255]/255,...
    [204 102 0]/255, [255 102 102]/255, [255 200 0]/255, [102 0 204]/255,...
    [153 0 76]/255, [0 153 76]/255, [20 20 255]/255,...
    [204 102 0]/255,...
    };
d = [1 2 2 2 2 2 2 2 1 1 1 1 2 2 2]; dashes = {'-', ':'};
lws = [3 1 1 1 1 1 1 1 2 2 2 2 1 1 1];
lws2 = [3 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
classN = {'Tin', 'Tout', 'Min'};

figure('Position', [600 100 800 900]); hold on
% Tin - Tout - Min
pltC = {[1 9 12 2 5 6 8],...
    [1 10 12 3 5 7 8],...
    [1 11 4 6 7 8],...
    };
pltC = {[1 9 12 2 5 6],...
    [1 10 12 3 5 7],...
    [1 11 4 6 7],...
    };
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*2+1); hold on
    legs = []; l = 0; 
    for sess = 1
        for e = pltC{c}
            plot(winCenter, smooth(pc_dec{e, sess}, sWin), 'Color', plotCols{e}, 'LineWidth', lws(e))
            l = l + 1; legs{l} = exclName{e};
        end
    end
    % legend(exclName)
    title([classN{c} ': ' par_sess{sess}.date])
    xlabel('time from motin onset')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    % legend(legs)
end
% -------------------------------
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*2+2); hold on
    legs = []; l = 0;
    for e = pltC{c}
        pltMat = nan(8, size(pc_dec{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_dec{e, sess}, 2)>0
            pltMat(sess, :) = pc_dec{e, sess};
            end
        end
        plot(winCenter, smooth(nanmean(pltMat, 1), sWin), 'Color', plotCols{e}, 'LineWidth', lws(e))
        l = l + 1; legs{l} = exclName{e};
    end
    % legend(exclName)
    title([classN{c} ': mean across sessions'])
    xlabel('time from motin onset')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs)
end
% ====================================
% Plot just's and no's separately
rows = 2; cols = 2; 
xlims = [0.1 0.5]; ylims = [0.45 0.75];
figure('Position', [600 100 800 900]); hold on
pltC = {[1 9 10 11 12],...
    [1 2 3 4 5],...
    };
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+1); hold on
    legs = []; l = 0;
    for sess = 1
        for e = pltC{c}
            plot(winCenter, smooth(pc_dec{e, sess}, sWin),...
                'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
            l = l + 1; legs{l} = exclName{e};
        end
    end
    % legend(exclName)
    title([par_sess{sess}.date])
    xlabel('time from motin onset')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end
% -------------------------------
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+2); hold on
    legs = []; l = 0;
    for e = pltC{c}
        pltMat = nan(8, size(pc_dec{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_dec{e, sess}, 2)>0
                pltMat(sess, :) = pc_dec{e, sess};
            end
        end
        plot(winCenter, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        l = l + 1; legs{l} = exclName{e};
    end
    title(['mean across sessions'])
    xlabel('time from motin onset')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end

%% Plot only just's and no's separately
% ====================================
% Plot just's and no's separately - same time point
rows = 2; cols = 2; 
figure('Position', [600 100 800 900]); hold on
xlims = [-0.5 0]; ylims = [0.5 1];
pltC = {[1 9 10 11 12],...
    [1 2 3 4 5],...
    };
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+1); hold on
    legs = []; l = 0;
    for sess = 1
        for e = pltC{c}
            plot(trl_dec, smooth(pc_decR{e, sess}, sWin),...
                'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
            l = l + 1; legs{l} = exclName{e};
        end
    end
    % legend(exclName)
    title([par_sess{sess}.date])
    xlabel('time to response')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end
% -------------------------------
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+2); hold on
    legs = []; l = 0;
    for e = pltC{c}
        pltMat = nan(8, size(pc_decR{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_decR{e, sess}, 2)>0
                pltMat(sess, :) = pc_decR{e, sess};
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        l = l + 1; legs{l} = exclName{e};
    end
    title(['mean across sessions'])
    xlabel('time to response')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end
suptitle('weights: same time point')
% ====================================
% Plot just's and no's separately - t=450ms stimulus-locked
rows = 2; cols = 2; 
figure('Position', [600 100 800 900]); hold on
xlims = [-0.5 0]; ylims = [0.5 1];
pltC = {[1 9 10 11 12],...
    [1 2 3 4 5],...
    };
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+1); hold on
    legs = []; l = 0;
    for sess = 1
        for e = pltC{c}
            plot(trl_dec, smooth(pc475_decR{e, sess}, sWin),...
                'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
            l = l + 1; legs{l} = exclName{e};
        end
    end
    % legend(exclName)
    title([par_sess{sess}.date])
    xlabel('time to response')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end
% -------------------------------
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+2); hold on
    legs = []; l = 0;
    for e = pltC{c}
        pltMat = nan(8, size(pc475_decR{1, 1}, 2));
        for sess = 1 : 8
            if size(pc475_decR{e, sess}, 2)>0
                pltMat(sess, :) = pc475_decR{e, sess};
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        l = l + 1; legs{l} = exclName{e};
    end
    title(['mean across sessions'])
    xlabel('time to response')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end
suptitle('weights: t=450ms')

% ====================================
% Plot just's and no's separately - same time point
rows = 2; cols = 2; 
figure('Position', [600 100 800 900]); hold on
xlims = [-0.5 0]; ylims = [0.5 1];
pltC = {[1 9 10 11 12],...
    [1 2 3 4 5],...
    };
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+1); hold on
    legs = []; l = 0;
    for sess = 1
        for e = pltC{c}
            plot(trl_dec, smooth(pc_decR{e, sess}, sWin),...
                'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
            l = l + 1; legs{l} = exclName{e};
        end
    end
    % legend(exclName)
    title([par_sess{sess}.date])
    xlabel('time to response')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end
% -------------------------------
for c = 1 : length(pltC)
    subplot(rows, cols, (c-1)*cols+2); hold on
    legs = []; l = 0;
    for e = pltC{c}
        pltMat = nan(8, size(pc_decR{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_decR{e, sess}, 2)>0
                pltMat(sess, :) = pc_decR{e, sess};
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        l = l + 1; legs{l} = exclName{e};
    end
    title(['mean across sessions'])
    xlabel('time to response')
    ylabel('xval decoder perf')
    xlim(xlims)
    ylim(ylims)
    legend(legs, 'Location', 'NorthWest')
    grid on
end
suptitle('weights: same time point')
% ====================================
% Plot just's and no's separately - SL & RL
rows = 2; cols = 3; 
figure('Position', [600 100 900 700]); hold on
xlims = [-0.5 0]; ylims = [0.45 1];
pltC = {[1 2 3 4 5 9 10 11 12],...
    };
xlimits = {[0 0.5], [-0.5 0], [-0.5 0],...
    [0 0.5], [-0.5 0], [-0.5 0]};
sess = 1; titles = {[par_sess{sess}.date ': SL'], [par_sess{sess}.date ': weights same TP'], [par_sess{sess}.date ': weights t=450'],...
    ['all sessions: SL'], ['all sessions: weights same TP'], ['all sessions: weights t=450']};
for c = 1 : length(pltC)
    legs = []; l = 0;
    for e = pltC{c}
        sess=1
        subplot(rows, cols, 1); hold on
        plot(winCenter, smooth(pc_dec{e, sess}, sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        subplot(rows, cols, 2); hold on
        plot(trl_dec, smooth(pc_decR{e, sess}, sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        subplot(rows, cols, 3); hold on
        plot(trl_dec, smooth(pc475_decR{e, sess}, sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        
        subplot(rows, cols, 4); hold on
        pltMat = nan(8, size(pc_dec{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_dec{e, sess}, 2)>0
                pltMat(sess, :) = pc_dec{e, sess};
            end
        end
        plot(winCenter, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        
        subplot(rows, cols, 5); hold on
        pltMat = nan(8, size(pc_decR{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_decR{e, sess}, 2)>0
                pltMat(sess, :) = pc_decR{e, sess};
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        
        subplot(rows, cols, 6); hold on
        pltMat = nan(8, size(pc475_decR{1, 1}, 2));
        for sess = 1 : 8
            if size(pc475_decR{e, sess}, 2)>0
                pltMat(sess, :) = pc475_decR{e, sess};
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        
        l = l + 1; legs{l} = exclName{e};
    end
    % legend(exclName)
    for plt = 1 : 6
        subplot(rows, cols, plt); hold on
        title(titles{plt})
        xlabel('time to response')
        ylabel('xval decoder perf')
        xlim(xlimits{plt})
        ylim(ylims)
        % legend(legs, 'Location', 'NorthWest')
        grid on
    end
end
legend(legs, 'Location', 'NorthWest')



% ====================================
% Plot just's and no's separately - SL & RL across session mean
rows = 1; cols = 2; 
figure('Position', [600 100 900 700]); hold on
xlims = [-0.5 0]; ylims = [0.45 1];
pltC = {[1 2 3 4 5 9 10 11 12],...
    };
xlimits = {[0 0.5], [-0.5 0], [-0.5 0],...
    [0 0.5], [-0.5 0], [-0.5 0]};
sess = 1; titles = {['all sessions: SL'], ['all sessions: weights same TP']};
for c = 1 : length(pltC)
    legs = []; l = 0;
    for e = pltC{c}
        subplot(rows, cols, 1); hold on
        pltMat = nan(8, size(pc_dec{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_dec{e, sess}, 2)>0
                pltMat(sess, :) = pc_dec{e, sess};
            end
        end
        plot(winCenter, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        xlabel('time after motion onset')
        
        subplot(rows, cols, 2); hold on
        pltMat = nan(8, size(pc_decR{1, 1}, 2));
        for sess = 1 : 8
            if size(pc_decR{e, sess}, 2)>0
                pltMat(sess, :) = pc_decR{e, sess};
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        xlabel('time to response')
        
        l = l + 1; legs{l} = exclName{e};
    end
    % legend(exclName)
    for plt = 1 : rows*cols
        subplot(rows, cols, plt); hold on
        title(titles{plt})
        ylabel('xval decoder perf')
        xlim(xlimits{plt})
        ylim(ylims)
        % legend(legs, 'Location', 'NorthWest')
        grid on
    end
end
legend(legs, 'Location', 'NorthWest')

%% Other plots
xlims = [-0.55 -0.02];
ylims = [0.45 1];
f = 4; pltFilt = ones(1, f)/f; wSmooth = 3;
figure; hold on
subplot(3, 2, 1); hold on
% e=1; plot(trl_dec, conv(pc_decR{e, sess}, pltFilt, 'same'), 'k', 'LineWidth', 2)
e=1; plot(trl_dec, smooth(pc_decR{e, sess}, wSmooth), 'k', 'LineWidth', 2)
for e = 2 : 8
%     plot(trl_dec, conv(pc_decR{e, sess}, pltFilt, 'same'))
    plot(trl_dec, smooth(pc_decR{e, sess}, wSmooth))
end
ylim(ylims); % xlim(xlims)
line(-0.1*[1 1], ylims, 'Color', 'k')
title('weights: same TP')

subplot(3, 2, 2); hold on
% e = 1; plot(trl_dec, conv(pc450_decR{e, sess}, pltFilt, 'same'), 'k', 'LineWidth', 2)
e = 1; plot(trl_dec, smooth(pc450_decR{e, sess}, wSmooth), 'k', 'LineWidth', 2)
for e = 2 : 8
%     plot(trl_dec, conv(pc450_decR{e, sess}, pltFilt, 'same'))
    plot(trl_dec, smooth(pc450_decR{e, sess}, wSmooth))
end
legend(exclName)
ylim(ylims)
line(-0.1*[1 1], ylims, 'Color', 'k')
title('weights: t=450')

for e = 1 : 8
    subplot(3, 4, e+4); hold on
    plot(trl_dec, smooth(pc_decR{e, sess}, wSmooth), 'LineWidth', 2)
    plot(trl_dec, smooth(pc450_decR{e, sess}, wSmooth))
%     plot(trl_dec, conv(pc_decR{e, sess}, pltFilt, 'same'), 'LineWidth', 2)
%     plot(trl_dec, conv(pc450_decR{e, sess}, pltFilt, 'same'))
    title(exclName{e})
    ylim(ylims)
    line(-0.1*[1 1], ylims, 'Color', 'k')
end
legend('weights: same TP', 'weights: t=450')
%% Plot RL decoding performance for all 8 sessions

doPlot = [1 1 0 0 1 1 0 1 1];
xlims = [-0.55 -0.02];
ylims = [0.45 1];
f = 4; pltFilt = ones(1, f)/f; wSmooth = 1;

figure; hold on
for sess = 1 : 8
    subplot(4, 2, sess); hold on
    legs = []; l = 0;
    % e=1; plot(trl_dec, conv(pc_decR{e, sess}, pltFilt, 'same'), 'k', 'LineWidth', 2)
    e=1; plot(trl_dec, smooth(pc_decR{e, sess}, wSmooth), 'k', 'LineWidth', 2)
    l = l + 1; legs{l} = exclName{e};
    for e = 2 : 8
        if doPlot(e) == 1
            %     plot(trl_dec, conv(pc_decR{e, sess}, pltFilt, 'same'))
            plot(trl_dec, smooth(pc_decR{e, sess}, wSmooth))
            l = l + 1; legs{l} = exclName{e};
        end
    end
    e=9; plot(trl_dec, smooth(pc_decR{e, sess}, wSmooth), ':r', 'LineWidth', 2)
    l = l + 1; legs{l} = exclName{e};
    ylim(ylims); % xlim(xlims)
    line(-0.1*[1 1], ylims, 'Color', 'k')
    title(par_sess{sess}.date)
    if sess > 6
        xlabel('time to response')
    end
    if rem(sess, 2) == 1
        ylabel('xval decoder perf')
    end
end
suptitle('weights: same TP')
legend(legs)

%---


% Example session and average of all sessions
wSmooth = 3;
figure; hold on
subplot(1, 2, 1); hold on
for sess = 1
    legs = []; l = 0;
    % e=1; plot(trl_dec, conv(pc_decR{e, sess}, pltFilt, 'same'), 'k', 'LineWidth', 2)
    e=1; plot(trl_dec, smooth(pc450_decR{e, sess}, wSmooth), 'k', 'LineWidth', 2)
    l = l + 1; legs{l} = exclName{e};
    for e = 2 : 8
        if doPlot(e) == 1
            %     plot(trl_dec, conv(pc_decR{e, sess}, pltFilt, 'same'))
            plot(trl_dec, smooth(pc450_decR{e, sess}, wSmooth))
            l = l + 1; legs{l} = exclName{e};
        end
    end
    e=9; plot(trl_dec, smooth(pc450_decR{e, sess}, wSmooth), ':r', 'LineWidth', 2)
    l = l + 1; legs{l} = exclName{e};
    ylim(ylims); % xlim(xlims)
    line(-0.1*[1 1], ylims, 'Color', 'k')
    title(par_sess{sess}.date)
    xlabel('time to response')
    ylabel('xval decoder perf')
end
subplot(1, 2, 2); hold on

legs = []; l = 0;
e=1;
pltMat = [];
for sess = 1 : 8
    pltMat(sess, :) = pc450_decR{e, sess};
end
plot(trl_dec, nanmean(pltMat, 1), 'k', 'LineWidth', 2)
l = l + 1; legs{l} = exclName{e};
for e = 2 : 8
    if doPlot(e) == 1
        pltMat = [];
        for sess = 1 : 8
            pltMat(sess, :) = pc450_decR{e, sess};
        end
        plot(trl_dec, nanmean(pltMat, 1))
        l = l + 1; legs{l} = exclName{e};
    end
end
e=9;
pltMat = [];
for sess = 1 : 8
    pltMat(sess, :) = pc450_decR{e, sess};
end
plot(trl_dec, nanmean(pltMat, 1), ':r', 'LineWidth', 2)
l = l + 1; legs{l} = exclName{e};
ylim(ylims); % xlim(xlims)
line(-0.1*[1 1], ylims, 'Color', 'k')
title('mean across sessions')
xlabel('time to response')
ylabel('xval decoder perf')
suptitle('weights: t=450ms')
legend(legs)
% =============
%% Use decoder derived at -100ms RL for all other time points

for sess = 1 : 8
    par = par_sess{sess};
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = par.beh.cho_trg;
    
    allNs = 1 : size(sl_sessMatZ{sess}, 1);
    a = ismember(allNs, par.unitIdxLIP_Tin);
    notTin = allNs(a == 0);
    
    excl = {[],...
        par.unitIdxLIP_Tin,...
        par.unitIdxLIP_Tout,...
        par.unitIdx_DinRFcC,...
        [par.unitIdxLIP_Tin par.unitIdxLIP_Tout],...
        [par.unitIdxLIP_Tin par.unitIdx_DinRFcC],...
        [par.unitIdxLIP_Tout par.unitIdx_DinRFcC],...
        [par.unitIdxLIP_Tin par.unitIdxLIP_Tout par.unitIdx_DinRFcC],...
        notTin,...
        };
    
    for e = 1 : length(excl)
        
        disp(num2str(sess))
        
        clear pcDD_test pDDtest_se useN
        
        useN = ones(size(rl_sessMatZ_short{sess}, 1), 1); useN(excl{e}) = 0; useN(sessZ{sess}.s==0) = 0;
        useN_decR{e, sess} = find(useN==1);
        
        w_proj = w_dec{e, sess}(i_proj, :);
        
        for i = 1:length(trl_dec)
            
            dW = find(trl_short >= trl_dec(i)-winSize/2 & trl_short <= trl_dec(i)+winSize/2); % Decoding window relative to dots onset
            
            xDD = squeeze(nanmean(rl_sessMatZ_short{sess}(useN==1, l, dW), 3));
            
            xDD = xDD - nanmean(xDD,2) ;
            lDD = (nanmean(~isnan(xDD), 1) == 1)';
            
            train = ~logical(mod(1:nL,2))' ;
            xD = xDD(:, train&lDD)' ;
            xTest = xDD(:, ~train&lDD)' ;
            c = yChoice(train&lDD);
            
            [wDD,sDD]     =  lassoglm(-xD, c, 'binomial','link','logit','Lambda',0.03);
            wDD = [0 ; wDD] ;

            choicePredTrain  = glmval(wDD,  -xD, 'logit') > 0.5;
            
            pcDD_train  =  sum(choicePredTrain==yChoice(train&lDD)) ./ sum(train&lDD) ;
            
            choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
            
            pcDD_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
            
            pDDtest_se(i) = sqrt( (pcDD_test(i)*(1-pcDD_test(i))) / sum(~train&lDD)) ;
            
            w_decR{e, sess}(i, :) = wDD;
            
            choicePredTest  = glmval(w_proj',  -xTest, 'logit') > 0.5;
            
            pcDD550_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
            
            pDDtest550_se(i) = sqrt( (pcDD550_test(i)*(1-pcDD550_test(i))) / sum(~train&lDD)) ;
            
        end
        pc_decR{e, sess} = pcDD_test;
        se_decR{e, sess} = pDDtest_se;
        
        pc475_decR{e, sess} = pcDD550_test;
        se475_decR{e, sess} = pDDtest550_se;
        
    end
end

%% =============================
%% Compare weights in different dimensions


for sess = 1 : 8
    useN = ones(size(rl_sessMatZ_short{sess}, 1), 1);
    allN = useN_dec{1, sess};
    
    for e1 = 1 : 12
        n1 = useN_dec{e1, sess};
        a1 = ismember(allN, n1);
        
        for e2 = 1 : 12
            n2 = useN_dec{e2, sess};
            a2 = ismember(allN, n2);
            overlap = allN(a1 & a2);
            
            if e1 ~= e2 & length(overlap) > 0
                
                if length(overlap) > 0
                    for i = 1 : size(w_dec{e1, sess}, 1)
                        
                        ww = w_dec{e1, sess}(i, 2:end);
                        nn = useN_dec{e1, sess}';
                        w1 = zeros(1, length(useN));
                        w1(nn) = ww; w1 = w1(overlap);
                        
                        ww = w_dec{e2, sess}(i, 2:end);
                        nn = useN_dec{e2, sess}';
                        w2 = zeros(1, length(useN));
                        w2(nn) = ww; w2 = w2(overlap);
                        
                        cs = acos(sum(w1.*w2) / (sqrt(sum(w1.^2))*sqrt(sum(w2.^2))));
                        % cs = acos(dot(w1, w2) / (sqrt(sum(w1.^2))*sqrt(sum(w2.^2))));
                        cs_dec(sess, e1, e2, i) = cs;
                        cs = (sum(w1.*w2) / (sqrt(sum(w1.^2))*sqrt(sum(w2.^2))));
                        cosim_dec(sess, e1, e2, i) = cs;
                        
                    end
                    
                    for i = 1 : size(w_decR{e1, sess}, 1)
                        
                        ww = w_decR{e1, sess}(i, 2:end);
                        nn = useN_dec{e1, sess}';
                        w1 = zeros(1, length(useN));
                        w1(nn) = ww;  w1 = w1(overlap);
                        
                        ww = w_decR{e2, sess}(i, 2:end);
                        nn = useN_dec{e2, sess}';
                        w2 = zeros(1, length(useN));
                        w2(nn) = ww;  w2 = w2(overlap);
                        
                        cs = acos(sum(w1.*w2) / (sqrt(sum(w1.^2))*sqrt(sum(w2.^2))));
                        % cs = acos(dot(w1, w2) / (sqrt(sum(w1.^2))*sqrt(sum(w2.^2))));
                        cs_decR(sess, e1, e2, i) = cs;
                        cs = (sum(w1.*w2) / (sqrt(sum(w1.^2))*sqrt(sum(w2.^2))));
                        cosim_decR(sess, e1, e2, i) = cs;
                        
                    end
                end
            end
        end
    end
end

% pdist([w1;w2], 'Cosine')

i_proj = find(round(winCenter, 3) == 0.476);

figure; hold on
imagesc(squeeze(cs_dec(1, :, :, i_proj)))
axis tight


figure; hold on
e1 = 1; e2 = 2;
for sess = 1 : 7
    plot(winCenter, smooth(real(squeeze(cosim_decR(sess, e1, e2, :))), 5))
end

% -------------

rows = 2; cols = 2; e1=1; sWin = 11;
figure('Position', [600 100 900 700]); hold on
ylims = [0 1];
pltC = {[2 3 4 5 9 10 11 12],...
    };
xlimits = {[0 0.5], [-0.5 0],...
    [0 0.5], [-0.5 0]};
ylims = [0.5 1];
sess = 1; titles = {[par_sess{sess}.date ': SL'], [par_sess{sess}.date ': weights same TP'],...
    ['all sessions: SL'], ['all sessions: weights same TP']};
for c = 1 : length(pltC)
    legs = []; l = 0;
    for e = pltC{c}
        sess=1; ylims = [0.5 1];
        subplot(rows, cols, 1); hold on
        plot(winCenter, smooth(real(squeeze(cosim_dec(sess, e1, e, :))), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        ylim(ylims)
        subplot(rows, cols, 2); hold on
        plot(trl_dec, smooth(real(squeeze(cosim_decR(sess, e1, e, :))), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        ylim(ylims)
        
        ylims = [0.6 1];
        subplot(rows, cols, 3); hold on
        pltMat = nan(8, size(cosim_dec, 4));
        for sess = 1 : 8
            if size(pc_dec{e, sess}, 2)>0
                pltMat(sess, :) = real(squeeze(cosim_dec(sess, e1, e, :)));
            end
        end
        plot(winCenter, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        ylim(ylims)
        
        subplot(rows, cols, 4); hold on
        pltMat = nan(8, size(cosim_decR, 4));
        for sess = 1 : 8
            if size(pc_decR{e, sess}, 2)>0
                pltMat(sess, :) = real(squeeze(cosim_decR(sess, e1, e, :)));
            end
        end
        plot(trl_dec, smooth(nanmean(pltMat, 1), sWin),...
            'Color', plotCols2{e}, 'LineWidth', lws2(e), 'LineStyle', dashes{d(e)})
        ylim(ylims)
        
        l = l + 1; legs{l} = exclName{e};
    end
    % legend(exclName)
    for plt = 1 : 4
        subplot(rows, cols, plt); hold on
        title(titles{plt})
        if rem(plt, 2) == 0; xlabel('time to response'); else
            xlabel('time after motion onset'); end
        ylabel('xval decoder perf')
        xlim(xlimits{plt})
        % legend(legs, 'Location', 'NorthWest')
        grid on
    end
end
legend(legs, 'Location', 'NorthWest')
suptitle('Cos Sim to decoder with all neurons')


save('E:\Matalb analyses\decoder_excluded',...
    'w_dec', 'pc_dec', 'se_dec',...
    'w_decR', 'pc_decR', 'se_decR',...
    'pc475_decR', 'se475_decR',...
    'cs_dec', 'cosim_dec', 'cs_decR', 'cosim_decR',...
    'excl', 'exclName', 'winCenter', 'trl_dec', 'i_proj')






