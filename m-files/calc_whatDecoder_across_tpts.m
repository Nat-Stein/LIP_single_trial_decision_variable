%% Across time-point accuracy

load(fullfile(saveLoc, 'par_sess'))
load(fullfile(saveLoc, 'decoder_excl_SL'), 'excl', 'exclName', 'useN_dec', 'winCenter', 'winStart', 'winEnd')

for sess = 1 : 8
    
    disp(num2str(sess))
    
    par = par_sess{sess};
    
    %% Load normalized activity for this session
    clear spkZ
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZ')
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt > 0.6; % l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = par.beh.cho_trg(1,l);
    
    % useN = ones(length(allNs), 1); useN(excl{e}) = 0; useN(spkZ.s==0) = 0;
    useNs = useN_dec{e, sess};
    
    clear pc_all_matrix se_all_matrix
    if length(useNs) > 0
        for i = 1 : length(winStart)
            
            dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
            
            xDD = squeeze(nanmean(spkZ.SL(useNs, l, dW), 3));
            xDD = xDD - nanmean(xDD,2) ;
            lDD = (nanmean(~isnan(xDD), 1) == 1)';
            
            train = ~logical(mod(1:nL,2))' ;
            xD = xDD(:, train&lDD)' ;
            c = yChoice(train&lDD);
            
            % Compute at time point
            [wDD,sDD]     =  lassoglm(-xD, c', 'binomial','link','logit','Lambda',0.03);
            wDD = [0 ; wDD] ;
            choicePredTrain  = glmval(wDD,  -xD, 'logit') > 0.5;
            pcDD_train  =  sum(choicePredTrain==yChoice(train&lDD)) ./ sum(train&lDD) ;
            
            % test at other time points
            for k = 1 : length(winStart)
                
                dW_k = find(par.tsl >= winStart(k) & par.tsl <= winEnd(k)); % Decoding window relative to dots onset
                
                xDD_k = squeeze(nanmean(spkZ.SL(useNs, l, dW_k), 3));
                xDD_k = xDD_k - nanmean(xDD_k,2) ;
                
                xTest = xDD_k(:, ~train&lDD)' ;
                c_k = yChoice(~train&lDD);
                
                choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
                
                pc_all_matrix(i, k)  =  sum(choicePredTest'==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                se_all_matrix(i, k) = sqrt( (pc_all_matrix(i, k)*(1-pc_all_matrix(i, k))) / sum(~train&lDD)) ;
                
            end
        end
    end
    pc_allMat(sess, :, :) = pc_all_matrix;
    
end
save(fullfile(saveLoc, 'decoderPerf_crossTP'), 'pc_allMat', 'winStart', 'winCenter', 'winEnd')
