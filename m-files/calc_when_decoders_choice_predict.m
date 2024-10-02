%% Computing When-decoders choice prediction
disp('Computing When decoder performance...')

load(fullfile(saveLoc, 'decoder_excl_SL'))

d = 6; % WhenD

clear pc_whenC_sign_dec se_whenC_sign_dec
for sess = 1 : 8
    disp(['Session ' num2str(sess) ': Computing When decoder performance...'])
    
    load(fullfile(saveLoc, [dNames{d} '_dat_S' num2str(sess)]),...
        'sl_rdm', 'rl_rdm', 'par')
    
    % par = par_sess{sess};
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt > 0.6;
    nL = sum(l);
    yChoice = par.beh.cho_trg(l);
    
    
    for i = 1 : length(winStart)
        
        dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
        
        % When-decoder sign test
        xDD_dim = nanmean(sl_rdm(l, dW), 2); xDD_dim = xDD_dim - nanmean(xDD_dim); lDD = (~isnan(xDD_dim) == 1)';
        cho_sign = sign(xDD_dim);
        choicePredSign = cho_sign(lDD) == -1;
        pcR_whenC_sign_test(i)  =  sum(choicePredSign==yChoice(lDD)') ./ sum(lDD) ;
        
        seR_whenC_sign_test(i) = sqrt( (pcR_whenC_sign_test(i)*(1-pcR_whenC_sign_test(i))) / sum(lDD)) ;
        
    end
    pc_whenC_sign_dec{sess} = pcR_whenC_sign_test;
    se_whenC_sign_dec{sess} = seR_whenC_sign_test;
    
    
end

disp('Saving stimulus-locked decoders...')
save(fullfile(saveLoc, 'decoder_excl_SL'), 'w_dec', 'pc_dec',...
    'se_dec', 'pc_whenC_sign_dec', 'se_whenC_sign_dec',...
    'excl_sess', 'exclName', 'useN_dec', 'winCenter', 'winStart', 'winEnd', 'winSize')
