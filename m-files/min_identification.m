%% Identify motion sensitive neurons based on passive motion viewing task

thFRvis = 1;
load(fullfile(saveLoc, ['spk_S' num2str(1)]), 'par')


fi = 120;  filtView = ones(1, fi)/fi;
t_vis = find(par.tsl >= 0 & par.tsl <= 0.08);
t_vis2 = find(par.tsl >= 0.05 & par.tsl <= 0.15);
t_ev = find(par.tsl >= 0.15 & par.tsl <= 1);
t_ev2 = find(par.tsl >= 0.3 & par.tsl <= 1);
t_bl = find(par.tsl >= -0.2 & par.tsl <= 0);
thAUC = 0.6; thFR = 2; 
doPlot = 0;

popView_SLnan_sess = [];
for sess = 1 : 8 
    
    clear popView_SLnan par    
     load(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkViewZnan', 'par')
    
    % load(fullfile(par.matPath, 'popView_SL'), 'popView_SLnan')
    popView_SLnan_sess{sess} = spkViewZnan.SL; % popView_SLnan;
    cohs = nanunique(par.behView.sig_coh);
    
    isMinI = nan(1, size(spkViewZnan.SL, 1)); isMinC = isMinI;
    isVin = isMinI; allAUC = isMinI; minCat = isMinI;
    pltCells = 1 : length (isMinI); % 
    rows = 3; colms = ceil(length(pltCells)/rows);
    if doPlot == 1; figure; hold on; plt = 0; end
    clear allAUC allAUC2 allBL allRamp allRampBL
    for n = pltCells
        pltFR = conv(squeeze(nanmean(spkViewZnan.SL(n, :, :), 2))'*1000, filtView, 'same');
        dFR = conv(diff(pltFR), ones(1, 10)/10, 'same');
        ramp = nanmean(dFR(1, t_vis));
        rampBL = nanmean(dFR(1, t_bl));
        ampBL = squeeze(nanmean(nanmean(spkViewZnan.SL(n, :, t_bl), 3), 2));
        ampVis = squeeze(nanmean(nanmean(spkViewZnan.SL(n, :, t_vis2), 3), 2));
        for c = 1 : 2
            tt{c} = find(par.behView.sig_coh == cohs(c));
            evFR(c) = squeeze(nanmean(nanmean(spkViewZnan.SL(n, tt{c}, t_ev), 3), 2));
        end
        motD = evFR(1) - evFR(2);
        
        labels = (par.behView.sig_coh > 0);
        motFR = squeeze(nanmean(spkViewZnan.SL(n, :, t_ev), 3));
        motFR2 = squeeze(nanmean(spkViewZnan.SL(n, :, t_ev2), 3));
        posclass = 1;
        [X,Y, T, AUC] = perfcurve(labels,motFR,posclass);
        [X2,Y2, T2, AUC2] = perfcurve(labels,motFR2,posclass);
        if doPlot == 1
            plt = plt + 1; subplot(rows, colms, plt); hold on;
            plot(X, Y); plot(X2, Y2);
            title([num2str(n) ': auc=' num2str(AUC) '/'  num2str(AUC2)])
        end
        
        allAUC(n) = AUC;
        allAUC2(n) = AUC2;
        allBL(n) = ampBL;
        allRamp(n) = ramp;
        allRampBL(n) = rampBL;
        if ramp > 0 & ramp > rampBL/1 & (nanmean(ampVis) > ampBL + 5*thFRvis/1000 | nanmean(ampVis) > ampBL * 5)
            if (AUC > thAUC & AUC2 > 1-thAUC) | AUC2 > thAUC & nanmean(motFR) > thFR/1000
                isMinI(n) = 1;
                minCat(n) = -1;
            elseif (AUC < 1-thAUC & AUC2 < thAUC) | AUC2 < 1-thAUC & nanmean(motFR) > thFR/1000
                isMinC(n) = 1;
                minCat(n) = 1;
            else
                isVin(n) = 1;
                minCat(n) = 0;
            end
        end
    end
    minCat_sess{sess} = minCat;
    
    mc_sess{sess}.allAUC = allAUC;
    mc_sess{sess}.allAUC2 = allAUC2;
    mc_sess{sess}.allBL = allBL;
    mc_sess{sess}.allRamp = allRamp;
    mc_sess{sess}.allRampBL = allRampBL;
    mc_sess{sess}.t_vis = t_vis;
    mc_sess{sess}.t_ev = t_ev;
    mc_sess{sess}.t_ev2 = t_ev2;
    mc_sess{sess}.thAUC = thAUC;
    mc_sess{sess}.filtView = filtView;
    
    
end

reverts = [];
non_reverts = [];
for sess = 1 : 8
    reverts(1, sess) = length(find(mc_sess{sess}.allAUC > 0.6 & mc_sess{sess}.allAUC2 < 0.4));
    reverts(2, sess) = length(find(mc_sess{sess}.allAUC < 0.4 & mc_sess{sess}.allAUC2 > 0.6));
    
    non_reverts(1, sess) = length(find(mc_sess{sess}.allAUC > 0.6));
    non_reverts(2, sess) = length(find(mc_sess{sess}.allAUC2 > 0.6));
    non_reverts(3, sess) = length(find(mc_sess{sess}.allAUC < 0.4));
    non_reverts(4, sess) = length(find(mc_sess{sess}.allAUC2 < 0.4));
end

% -------------------
% Visual validation
plotMin = 0;
if plotMin == 1
    for sess = 1 : 8
        xlims = [-0.1 0.5];
        minCat = minCat_sess{sess};
        
        pltCells = find(isnan(minCat));
        rows = 5; colms = ceil(length(pltCells)/rows);
        par = par_sess{sess};
        cohs = unique(par.behView.sig_coh);
        
        figure('Position', [100 100 800 500]); hold on;
        for c = 1 : length(cohs)
            
            plt = 0;
            for n = pltCells
                plt = plt + 1;
                
                trls = find(par.behView.sig_coh == cohs(c));
                medRT = nanmedian(par.dataView(:, par.idx.dot_off)-par.dataView(:, par.idx.dot_on));
                pltTPT = find(par.tsl <= (medRT  - fw/2/1000));
                
                subplot(rows, colms, plt); hold on;
                title([num2str(n) ': ' num2str(minCat(n)) ', auc=' num2str(round(allAUC(n), 2)) '/'  num2str(round(allAUC2(n), 2))])
                pltFR = conv(squeeze(nanmean(popView_SLnan_sess{sess}(n, trls, :), 2))'*1000, filtView, 'same');
                plot(par.tsl(pltTPT), pltFR(pltTPT), 'Color', collies{c})
                
                % ------
                % Decision
                trls = find(par.dataMat(:, par.idx.sig_coh) == cohs(c));
                medRT = nanmedian(par.dataMat(trls, par.idx.rt))/1000;
                pltTPT = find(par.tsl <= (medRT - 0.1));
                
                subplot(rows, colms, plt); hold on;
                pltFR = conv(squeeze(nanmean(sl_sessMat0{sess}(n, trls, :), 2))'*1000, filtView, 'same');
                plot(par.tsl(pltTPT), pltFR(pltTPT),...
                    '--', 'Color', collies{c}, 'LineWidth', lw2)
            end
        end
        for plt = 1 : rows*colms
            subplot(rows, colms, plt); hold on;
            xlim(xlims)
            xlabel('Time after motion onset [s]')
            if plt == 1 | plt == 4
                ylabel('FR [sp/s]')
            else
                ylabel('[a. u.]')
            end
            set(gca, 'TickDir', 'Out', 'FontSize', 8)
        end
        suptitle(par.date)
    end
end
%% -------------------
%% Recategorize neurons based on visual inspection

Min_sess = [];
for sess = 1 : 8
    %par_sess{sess}.unitIdx_DinRFcC = par_sess{sess}.unitIdx_DinRFc;
    %par_sess{sess}.unitIdx_DinRFiC = par_sess{sess}.unitIdx_DinRFi;
    %par_sess{sess}.unitIdx_DinRFnC = par_sess{sess}.unitIdx_DinRFn;
    minCat_sessC{sess} = minCat_sess{sess};
end
sess=1;
minCat_sessC{sess}([]) = 1;
minCat_sessC{sess}([]) = -1;
minCat_sessC{sess}([77]) = 0;
minCat_sessC{sess}([53 62 174 119]) = NaN;
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);
sess= 2;
minCat_sessC{sess}([86 11 65]) = 1;
minCat_sessC{sess}([3]) = -1;
minCat_sessC{sess}([41 39]) = 0;
minCat_sessC{sess}([]) = NaN;
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);
sess=3;
minCat_sessC{sess}([]) = 1;
minCat_sessC{sess}([]) = -1;
minCat_sessC{sess}([30]) = 0;
minCat_sessC{sess}([9]) = NaN;
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);
sess=4;
minCat_sessC{sess}([4 97]) = 1;
minCat_sessC{sess}([95]) = -1;
minCat_sessC{sess}([45 51 114 132]) = 0; 
minCat_sessC{sess}([]) = NaN;
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);
sess=5;
minCat_sessC{sess}([27]) = 1;
minCat_sessC{sess}([63 50 93]) = -1;
minCat_sessC{sess}([43 65]) = 0;
minCat_sessC{sess}([]) = NaN;
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);
sess=6;
minCat_sessC{sess}([]) = 1;
minCat_sessC{sess}([104]) = -1;
minCat_sessC{sess}([9 133 135]) = 0;
minCat_sessC{sess}([7 20 100 115 28 ]) = NaN;
Min_sess{sess}.unitIdx_DinRFcNegC = [40];
Min_sess{sess}.unitIdx_DinRFiNegC = [38];
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);
sess=7;
minCat_sessC{sess}([16 27 55 101 106 108]) = 1;
minCat_sessC{sess}([6 140 ]) = -1;
minCat_sessC{sess}([23]) = 0;
minCat_sessC{sess}([8 123 129]) = NaN;
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);
sess=8;
minCat_sessC{sess}([112 8 54 127 174 186]) = 1;
minCat_sessC{sess}([115 116 143 152 29 32 146]) = -1;
minCat_sessC{sess}([159 51 82 102 144 183 188 196]) = 0; %199
minCat_sessC{sess}([36 55 56 65 70]) = NaN;
Min_sess{sess}.unit_MinC = find(minCat_sessC{sess} == 1);
Min_sess{sess}.unit_MinI = find(minCat_sessC{sess} == -1);
Min_sess{sess}.unit_MinN = find(minCat_sessC{sess} == 0);

for sess = 1 : 8
    Min_sess{sess}.minCat_sessC = minCat_sessC{sess};
    Min_sess{sess}.minCat_sess = minCat_sess{sess};
end

parMin.t_vis = t_vis;
parMin.t_vis2 = t_vis2;
parMin.t_ev = t_ev;
parMin.t_ev2 = t_ev2;
parMin.t_bl = t_bl;
parMin.thAUC = thAUC;
parMin.thFR = thFR;
parMin.filtView = filtView;
%% Remove neurons that are also Tin or Tout


load(fullfile(saveLoc, 'par_sess'))

for sess = 1 : 8
    
    par = par_sess{sess};
    
    minc = Min_sess{sess}.unit_MinC;
    p=ismember(minc, [par.unit_Tin par.unit_Tout]);
    minc=minc(~p);
    par_sess{sess}.unit_MinC = minc;
    
    minc =  Min_sess{sess}.unit_MinI;
    p=ismember(minc, [par.unit_Tin par.unit_Tout]);
    minc=minc(~p);
    par_sess{sess}.unit_MinI = minc;
    
    minc = Min_sess{sess}.unit_MinN;
    p=ismember(minc, [par.unit_Tin par.unit_Tout]);
    minc=minc(~p);
    par_sess{sess}.unit_MinN = minc;
    
    w = zeros(1, par.numCells);
    w(par_sess{sess}.unit_MinC) = 1/length(par_sess{sess}.unit_MinC); 
    if size(w, 1) > size(w, 2); w = w'; save(ff, 'w'); end
    ff = fullfile(saveLoc,...
        ['MinC_w_S' num2str(sess)]);
    save(ff, 'w')

    w = zeros(1, par.numCells);
    w(par_sess{sess}.unit_MinI) = 1/length(par_sess{sess}.unit_MinI); 
    if size(w, 1) > size(w, 2); w = w'; save(ff, 'w'); end
    ff = fullfile(saveLoc,...
        ['MinI_w_S' num2str(sess)]);
    save(ff, 'w')
    
end

save(fullfile(saveLoc, 'par_sess'), 'par_sess')

% clear overlap overlapCell
% for sess = 1 : 8
%     par = par_sess{sess};
%     sigs{1} = par.unit_Tin;
%     sigs{2} = par.unit_MinC;
%     sigs{3} = par.unit_Tout;
%     sigs{4} = par.unit_MinI;
%     
%     for s1 = 1 : 4
%         for s2 = 1 : 4
%             if s1 ~= s2
%                 p=ismember(sigs{s1}, sigs{s2});
%                 overlap(sess, s1, s2) = sum(p);
%                 overlapCell{sess, s1, s2} = sigs{s1}(p);
%                 
%             end
%         end
%     end
% end
% cname = {'TinC', 'MinC', 'TinI', 'MinI'};
% for sess = 1 : 8
%     disp('---------------------------------')
%    disp(num2str(sess)) 
%    for s1 = 1 : 4
%         for s2 = 1 : 4
%             if s1 < s2 && length(overlapCell{sess, s1, s2}) > 0
%                 
%                 for cc = overlapCell{sess, s1, s2}
%                    disp(['Cell ' num2str(cc) ': ' cname{s1} ' & ' cname{s2}]) 
%                 end
%                 
%             end
%         end
%    end
% end
% 
% newNeuronClasses = [];
% for sess = 1 : 8
%     par = par_sess{sess};
%     newNeuronClasses{sess}.unit_Tin = par.unit_Tin;
%     newNeuronClasses{sess}.unit_Tout = par.unit_Tout;
%     newMinC = [];
%     p=ismember(par.unit_MinC, [par.unit_Tout par.unit_Tin]);
%     newMinC = par.unit_MinC(p==0);
%     newNeuronClasses{sess}.unit_MinC = newMinC;
%     newNeuronClasses{sess}.unit_MinI = par.unit_MinI;
%     par_sess{sess}.unit_MinC = newNeuronClasses{sess}.unit_MinC;
%     newNeuronClasses{sess}.date = par.date;
% end
% save('E:\Matalb analyses\newNeuronClasses', 'newNeuronClasses')


        


