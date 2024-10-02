% Renamed from calc_Fig4a_S9


%% Decoder performace: fixed time window
%% weights from 450ms applied to all time points for 
%% 1. all neurons, 
%% 2. only Tin neurons, 
%% 3. no Tin neurons, 
%% 4. no Tin or Min neurons



load(fullfile(saveLoc, 'par_sess'))

%% Aligned to motion onset

% Take out the exclusion examples we're not using

load(fullfile(par_sess{1}.saveLoc, 'decoder_excl_SL'))
i_proj = find(round(winCenter, 3) == 0.451);

for sess = 1 : 8
    par = par_sess{sess};
    
    % Load normalized activity for this session
    clear spkZ
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZ')
    
    rt = par.beh.rt;
    coh = par.beh.dot_coh;
    l = rt > 0.6; % l = rt>0.6  & coh >= 0;
    nL = sum(l);
    yChoice = par.beh.cho_trg(l);
    
    for e = 1 : length(excl_sess{sess})
        
        disp(['Sess ' num2str(sess) 'e ' num2str(e) '/' num2str(length(excl))])
        
        clear pcDD_test pDDtest_se useN pDDtest450_se pcDD450_test
        
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
            
        end
    end
end

save(fullfile(par_sess{1}.saveLoc, 'decoder_excluded_450'), 'w_dec', 'pc_dec', 'se_dec', 'excl', 'exclName', 'winCenter',...
    'pc450_dec', 'se450_dec')










