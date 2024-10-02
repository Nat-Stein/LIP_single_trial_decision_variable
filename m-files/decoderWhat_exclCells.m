% LIP choice decoder including/excluding Tin, Tout, Min populations

% Use this one

% load('E:\Matalb analyses\popMatZ_sessSL1-10', 'sl_sessMatZ', 'sessZ')

% winSize = 0.001 ; 
% winStart = .0005:.005:.6-winSize;
% winEnd = winStart+winSize ;
% winCenter = mean([winStart ; winEnd]);

% rng default
% X = randn(100,10);
% weights = [0.6;0.5;0.7;0.4];
% y = X(:,[2 4 5 7])*weights + randn(100,1)*0.1; % Small added noise
% 
% B = lassoglm(X,y);

winSize = 0.05 ; 
winCenter = 0.001 : 0.005 : 0.6;
winStart = winCenter-winSize/2;
winEnd = winCenter+winSize/2;
winCenter = mean([winStart ; winEnd]);


exclName = {'All', 'noTin', 'noTout', 'noMinC', 'noTiTo', 'noTiMiC', 'noToMiC', 'noTiToMiC',...
    'justTin', 'justTout', 'justMinC', 'justTiTo', 'justTinMin', 'justMinCI', 'noTinMin', 'noMinCI'};

clear pc_dec se_dec w_dec

opts = statset('lassoglm');
opts.UseParallel = true;
opts.MaxIter = 1000;

load(fullfile(saveLoc, 'par_sess'))

for sess = 1 : 2
    par = par_sess{sess};
    
    load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZ', 'par')
    load(fullfile(saveLoc, ['beh_S' num2str(sess)]), 'beh')
    
    
    rt = beh.rt;
    coh = beh.dot_coh;
    l = rt>0.6; % l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = 1-beh.cho_trg;
    
    allNs = 1 : size(spkZ.SL, 1);
    a = ismember(allNs, par.unitIdxLIP_Tin);
    notTin = allNs(a == 0);
    a = ismember(allNs, par.unitIdxLIP_Tout);
    notTout = allNs(a == 0);
    a = ismember(allNs, par.unitIdx_DinRFcC);
    notMinC = allNs(a == 0);
    a = ismember(allNs, [par.unitIdxLIP_Tin par.unitIdxLIP_Tout]);
    notTinTout = allNs(a == 0);
    a = ismember(allNs, [par.unitIdxLIP_Tin par.unitIdxLIP_Tout par.unitIdx_DinRFcC par.unitIdx_DinRFiC]);
    notTinMin = allNs(a == 0);
    a = ismember(allNs, [par.unitIdx_DinRFcC par.unitIdx_DinRFiC]);
    notMinCI = allNs(a == 0);
    
    excl = {[],...
        par.unitIdxLIP_Tin,...
        par.unitIdxLIP_Tout,...
        par.unitIdx_DinRFcC,...
        [par.unitIdxLIP_Tin par.unitIdxLIP_Tout],...
        [par.unitIdxLIP_Tin par.unitIdx_DinRFcC],...
        [par.unitIdxLIP_Tout par.unitIdx_DinRFcC],...
        [par.unitIdxLIP_Tin par.unitIdxLIP_Tout par.unitIdx_DinRFcC],...
        notTin,...
        notTout,...
        notMinC,...
        notTinTout,...
        notTinMin,...
        notMinCI,...
        [par.unitIdxLIP_Tin par.unitIdxLIP_Tout par.unitIdx_DinRFcC par.unitIdx_DinRFiC],...
        [par.unitIdx_DinRFcC par.unitIdx_DinRFiC],...
        };
    
    excl_sess{sess} = excl;
    
    for e = 1 : length(excl)
        
        disp(['S' num2str(sess) ', ' num2str(e) '/' num2str(length(excl))])
        
        clear pcDD_test pDDtest_se useN
        
        useN = ones(size(spkZ.SL, 1), 1); useN(excl{e}) = 0; useN(spkZ.s==0) = 0;
        useN_dec{e, sess} = find(useN==1);
        
        if sum(useN) > 0
            for i = 1:length(winStart)
                
                dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
                
                xDD = squeeze(nanmean(spkZ.SL(useN==1, l, dW), 3));
                
                xDD = xDD - nanmean(xDD,2) ;
                lDD = (nanmean(~isnan(xDD), 1) == 1)';
                
                train = ~logical(mod(1:nL,2))' ;
                xD = xDD(:, train&lDD)' ;
                xTest = xDD(:, ~train&lDD)' ;
                c = yChoice(train&lDD); 
                
%                 c = [ones(size(c));c]; % c = [c; ones(size(c))];
                [wDD,sDD]     =  lassoglm(-xD, c', 'binomial','link','logit','Lambda',0.03);
                
                wDD = [0 ; wDD] ;
                
                choicePredTrain  = glmval(wDD,  -xD, 'logit') > 0.5;
                
                pcDD_train  =  sum(choicePredTrain==yChoice(train&lDD)) ./ sum(train&lDD) ;
                
                choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
                
                pcDD_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)') ./ sum(~train&lDD) ;
                
                pDDtest_se(i) = sqrt( (pcDD_test(i)*(1-pcDD_test(i))) / sum(~train&lDD)) ;
                
                w_dec{e, sess}(i, :) = wDD;
                
            end
            pc_dec{e, sess} = pcDD_test;
            se_dec{e, sess} = pDDtest_se;
        end
    end
end

figure; hold on
for e = 1 : length(excl)
    plot(winStart, pc_dec{e, sess})
end



save(saveLoc, 'w_dec', 'pc_dec', 'se_dec', 'excl', 'exclName')

%n3v3r20rt2p1k32




%% Grab weights at 450ms and make them into dimensions - stimulus aligned

i_proj = find(round(winCenter, 3) == 0.451); % i_proj = find(round(winCenter, 3) == 0.476); 
% i_proj = find(round(winCenter, 3) == 0.5510); 

for sess = 1 : 8
    clear useDat
    for tp = 1:length(winStart)
        tpts = find(par.tsl >= winStart(tp) & par.tsl <= winEnd(tp));
        useDat(:, :, tp) = nanmean(sl_sessMatZ{sess}(:, :, tpts), 3);
    end
    
    for e = 1 : size(pc_dec, 1)
        
        if size(w_dec{e, sess}, 1) >0
            w = w_dec{e, sess}(i_proj, 2:end); % used to be: w = w_dec{e, sess}(i_proj, :);
            %goodN = ones(size(rl_sessMatZ_short{sess}, 1), 1); goodN(excl{e}) = 0; goodN(sessZ{sess}.s==0) = 0;
            goodN = useN_dec{e, sess};
            
            C = squeeze(mmx('mult', w, useDat(goodN, :, :)));
            
            sl_exDec{e, sess} = C;
        end
    end
end

%% Decoder performace with weights from 450ms applied - stimulus locked


for sess = 1 : 8
    par = par_sess{sess};
    
    rt = beh.rt);
    coh = beh.dot_coh);
    l = rt>600; % l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = beh.cho_trg;
    
    for e = 1 : length(excl_sess{sess})
        
        disp(num2str(sess))
        
        
        clear pcDD_test pDDtest_se useN pDDtest550_se pcDD550_test
        
        useN = ones(size(rl_sessMatZ_short{sess}, 1), 1); useN(excl_sess{sess}{e}) = 0; useN(sessZ{sess}.s==0) = 0;
        
        if size(w_dec{e, sess}, 1) > 0
            
            w_proj = w_dec{e, sess}(i_proj, :);
            for i = 1:length(winStart)
                
                dW = find(par.tsl >= winStart(i) & par.tsl <= winEnd(i)); % Decoding window relative to dots onset
                
                
                xDD = squeeze(nanmean(sl_sessMatZ{sess}(useN==1, l, dW), 3));
                
                xDD = xDD - nanmean(xDD,2) ;
                lDD = (nanmean(~isnan(xDD), 1) == 1)';
                
                train = ~logical(mod(1:nL,2))' ;
                %xD = xDD(:, train&lDD)' ;
                xTest = xDD(:, ~train&lDD)' ;
                c = yChoice(train&lDD);
                
                choicePredTest  = glmval(w_proj',  -xTest, 'logit') > 0.5;
                
                pcDD550_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                pDDtest550_se(i) = sqrt( (pcDD550_test(i)*(1-pcDD550_test(i))) / sum(~train&lDD)) ;
                
            end
            pc450_dec{e, sess} = pcDD550_test;
            se450_dec{e, sess} = pDDtest550_se;
        end
    end
end


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
clear w_decR
for sess = 1 : 8
    par = par_sess{sess};
    
    rt = beh.rt);
    coh = beh.dot_coh);
    l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = beh.cho_trg;
    
%     allNs = 1 : size(rl_sessMatZ_short{sess}, 1);
%     a = ismember(allNs, par.unitIdxLIP_Tin);
%     notTin = allNs(a == 0);
%     a = ismember(allNs, par.unitIdxLIP_Tout);
%     notTout = allNs(a == 0);
%     a = ismember(allNs, par.unitIdx_DinRFcC);
%     notMin = allNs(a == 0);
%     a = ismember(allNs, [par.unitIdxLIP_Tin par.unitIdxLIP_Tout]);
%     notTinTout = allNs(a == 0);
%     
%     excl = {[],...
%         par.unitIdxLIP_Tin,...
%         par.unitIdxLIP_Tout,...
%         par.unitIdx_DinRFcC,...
%         [par.unitIdxLIP_Tin par.unitIdxLIP_Tout],...
%         [par.unitIdxLIP_Tin par.unitIdx_DinRFcC],...
%         [par.unitIdxLIP_Tout par.unitIdx_DinRFcC],...
%         [par.unitIdxLIP_Tin par.unitIdxLIP_Tout par.unitIdx_DinRFcC],...
%         notTin,...
%         notTout,...
%         notMin,...
%         notTinTout,...
%         };
    
    for e = 1 : length(excl)
        
        disp(num2str(sess))
        
        clear pcDD_test pDDtest_se useN
        
        useN = ones(size(rl_sessMatZ_short{sess}, 1), 1); useN(excl_sess{sess}{e}) = 0; useN(sessZ{sess}.s==0) = 0;
        useN_decR{e, sess} = find(useN==1);
        if sum(useN) > 0
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
                
                %[wDD,devDD,sDD]    =  glmfit(-xD, c, 'binomial','link','logit','Constant','Off');
                %[wChoice,sAll]     =  lassoglm(-xA, c, 'binomial','link','logit','Lambda',0.03);
                [wDD,sDD]     =  lassoglm(-xD, c, 'binomial','link','logit','Lambda',0.03);
                wDD = [0 ; wDD] ;
                
                %choicePredTrain  = glmval(wDD,  -xD, 'logit','Constant','Off') > 0.5;
                choicePredTrain  = glmval(wDD,  -xD, 'logit') > 0.5;
                
                pcDD_train  =  sum(choicePredTrain==yChoice(train&lDD)) ./ sum(train&lDD) ;
                
                % choicePredTest  = glmval(wDD,  -xTest, 'logit','Constant','Off') > 0.5;
                choicePredTest  = glmval(wDD,  -xTest, 'logit') > 0.5;
                
                pcDD_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                pDDtest_se(i) = sqrt( (pcDD_test(i)*(1-pcDD_test(i))) / sum(~train&lDD)) ;
                
                w_decR{e, sess}(i, :) = wDD;
                
                %choicePredTest  = glmval(w_proj',  -xTest, 'logit','Constant','Off') > 0.5;
                choicePredTest  = glmval(w_proj',  -xTest, 'logit') > 0.5;
                
                pcDD550_test(i)  =  sum(choicePredTest==yChoice(~train&lDD)) ./ sum(~train&lDD) ;
                
                pDDtest550_se(i) = sqrt( (pcDD550_test(i)*(1-pcDD550_test(i))) / sum(~train&lDD)) ;
                
            end
            pc_decR{e, sess} = pcDD_test;
            se_decR{e, sess} = pDDtest_se;
            
%             pc475_decR{e, sess} = pcDD550_test;
%             se475_decR{e, sess} = pDDtest550_se;
            
            pc450_decR{e, sess} = pcDD550_test;
            se450_decR{e, sess} = pDDtest550_se;
            
            %         pc550_decR{e, sess} = pcDD550_test;
            %         se550_decR{e, sess} = pDDtest550_se;
        end
    end
end


%% Plot stimulus- and response-aligned traces for the same time point
rows = 1; cols = 2; 
xLimsSL = [0 0.6]; 
xLimsRL = [-0.3 0]; 
ylims = [0.45 1];
pltC = {[1 9 14 12 13 2 16 5 15],... % [1 2 4 5 9 11 12 13],...
    [1 2 3 4 5],...
    };

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

[200 50 53]
[255 175 143]
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



%% Plot stimulus- and response-aligned traces for the same time point and weights derived at 450ms
rows = 1; cols = 2; 
xLimsSL = [0 0.6]; 
xLimsRL = [-0.3 0]; 
ylims = [0.45 1];
pltC = {[1],... % [1 2 5 9 11 12 13],...
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
legend('moving window', 'fixed window', 'Location', 'NorthWest')
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
%% OTHER PLOTS

%% Plot simulated data

aaa = load('C:\Users\shadlenlab\Dropbox\Natalie_and_Mike\sim4Natalie-29Apr')




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
    
    rt = beh.rt;
    coh = beh.dot_coh;
    l = rt>600  & coh >= 0;
    nL = sum(l);
    yChoice = beh.cho_trg;
    
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






