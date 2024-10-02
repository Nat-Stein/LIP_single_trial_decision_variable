% Plot stimulus-locked mean traces for all coherences
% meanTracesSL_Fig4


absCohs = unique(beh.dot_coh); absCohs = flip(absCohs);
dirs = unique(beh.dot_dir);
diros = [nanmean(beh.dot_dir(beh.dot_dir<180)),...
    nanmean(beh.dot_dir(beh.dot_dir>180))]; dirs = diros;
chos = unique(beh.cho_trg); if length(par.date) == 6; chos = chos(chos>0); end

if ~exist('plotLW'); plotLW = {0.5, 0.5}; end
if exist('changeCol'); colSet = changeCol; else; colSet{1} = plotColors_Coh; colSet{2} = plotColors_Coh; end


for coh = length(absCohs) : -1 : 1
    pltFR = nan(1, length(t_in));
    for d = 1 : 2
        if absCohs(coh) ~= 0
            trls = find(round(beh.dot_dir/180) == round(dirs(d)/180) &...
                beh.task == he &...
                ~isnan(beh.rt) &...
                beh.task == beh.effector &...
                beh.dot_coh == absCohs(coh));
        else
            trls = find(beh.task == he &...
                beh.task == beh.effector &...
                ~isnan(beh.rt) &...
                beh.dot_coh == absCohs(coh));
        end
        
        pltMat = (1-2*invDir) * squeeze(Cin(trls, :));
        if doMeanSub == 1
            pltMat = pltMat - repmat(squeeze(nanmean(Cin, 1)), length(trls), 1);
        end
        plotTPT = 1 : size(pltMat, 2);
        pltFR(1, plotTPT) = nanmean(pltMat(:, plotTPT), 1);
        convPLT = conv(squeeze(pltFR), plotFilt, 'same'); convPLT = convPLT(1,1:size(pltFR,2));
        if doBLcorr == 1
            bl(coh, d) = nanmean(convPLT(:, tBLv));
            convPLT = convPLT - bl(coh, d); %nanmean(nanmean(allDim_sl{dim}{sc, sess}(trls, tBLv), 2), 1)*1000;
        end
        thisCol = colSet{d}{coh};
        medianRT = median(beh.rt(trls))/1000;
        plotTPT = find(t_in <= medianRT - rtAdd);
        plot(t_in(plotTPT), convPLT(plotTPT),...
            'Color', thisCol, 'LineStyle', plotDash{3-d}, 'LineWidth', plotLW{3-d})
        maxTPT = [max(plotTPT(1), find(round(t_in, 3) == round(xLimsSL(1), 3))) :...
            min(plotTPT(end), find(round(t_in, 3) == round(xLimsSL(2), 3)))];
        maxY = max(maxY, max(convPLT(maxTPT)));
        minY = min(minY, min(convPLT(maxTPT)));
    end
end

xlim(xLimsSL)
xlabel('Time from motion onset [s]')
ylabel('Activity [a.u.]')
set(gca,'TickDir','out', 'FontSize', 10);
