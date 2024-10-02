% plot_singleTrials_highlights

%% Plot settings

%% Single trials Tin-c, PC1, Ramp
sess = 1; % Example trials taken from session 1
set(gcf, 'renderer', 'Painters');
% for sig = 1 : length(pltSig)
for d = 1 : length(plotSig)
    
    dim = plotSig{d};
    
    clear sl
    assignSignal
    clear rl
    
    maxY = 0; minY = 0;
    slMat = sl;
    if doMeanSub0 == 1
        trls = find(sigSL.sig_coh == 0);
        slMat = slMat - repmat(nanmean(slMat(trls, :), 1), size(slMat, 1), 1);
    end
    
    % Stimulus-locked
    plt = 0;
    for coh = plotCohs
        plt = plt + 1;
        subplot(rows, colms, (d-1)*colms + plt); hold on
        trls = find(sigSL.sig_coh == sigCohs(coh) &...
            sigSL.session == sess);
        k =  1.5 - sign(sigCohs(coh)) * 0.5;
        if sigCohs(coh) == 0
            thisCol = 0.1 * [1 1 1];
        else
            thisCol = colorsYeBl{round(k)}{coh};
        end
        
        pltMat = squeeze(slMat(trls, :));
        if meanCorr == 1
            pltMat = pltMat - repmat(nanmean(pltMat, 1), size(pltMat, 1), 1);
        end
        convPLT = conv2(1, plotFiltST, pltMat, 'same');
        plotTPT = 1 : length(convPLT);
        ct = 0;
        for trl = [1 : numTr uTRL{coh}]
            if length(find(exclTRL{coh} == trl)) == 0
                ct = ct + 1;
                
                plotVec = convPLT(trl, plotTPT);
                if doBLcorr == 1
                    plotVec = plotVec - repmat(nanmean(plotVec(:, tBLv), 2), 1, size(slMat, 2));
                    %yl = (yLims{d} - diff(yLims{d})/2);
                    yl = yLims_BLsub{d};
                else
                    yl = yLims_raw{d};
                end
                pltCol = thisCol;
                alpha = 0.2; lw = 0.5;
                if length(find(uTRL{coh} == trl)) > 0; alpha = 0.9; lw = 1.2; end
                pltCol = [pltCol alpha];
                plot(sigSL.t(plotTPT), plotVec,...
                    'Color', pltCol, 'LineStyle', pltDash{1}, 'LineWidth', lw)
                maxTPTS = find(sigSL.t >= xlims{1}(1) & sigSL.t <= xlims{1}(2));
                maxTPT = [max(plotTPT(1), find(round(sigSL.t, 3) == round(xlims{1}(1), 3))) :...
                    min(plotTPT(end), find(round(sigSL.t, 3) == round(xlims{1}(2), 3)))];
                maxY = max(maxY, max(max(convPLT(:, maxTPT))));
                minY = min(minY, min(min(convPLT(:, maxTPT))));
            end
        end
        
        text(0.22, yl(2) - 0.05 * diff(yl), cohText{coh}, 'FontSize', fontSize)
        
        xlim(xlims{1});
        ylim(yl)
        set(gca,'TickDir','out', 'FontSize', fontSize);
        
        yticks(yTicks{d})
        if coh == plotCohs(1)
            ylabel(sigName{d}, 'Fontsize', 9)
            hYlabel = get(gca, 'YLabel');
            set(hYlabel, 'rotation', 0, 'HorizontalAlignment','right')
        else
            ax = gca;
            numTicks = length(ax.YAxis.TickValues);
            for nt = 1 : numTicks; ax.YAxis.TickLabels{nt, 1} = ''; end
        end
        
        
        xticks([.2 .3 .4 .5]);
        if d == length(plotSig)
            xlabel('Time from motion onset [s]');
        else
            xticklabels({'', '','', ''})
            if d == 1
                yPos = yl(2) + diff(yl)*0.1;
                text(0.18, yPos, '\it Motion onset', 'FontSize', 9)
            end
        end
        
    end
end



