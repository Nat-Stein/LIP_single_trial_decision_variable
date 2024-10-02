% plot_trialAvg_decision_allSess
% new version of plot_bigSbigT_multipleS_trialAvg_decision_normS

singleSess = 1;
pSL = colms - plot_resp_aligned;
i = 0;

for d = 1 : length(plotSig)
    
    dim = plotSig{d};
    i = i + 1;
    
    assignSignal
    
    % Overwrite for single session plots
    ylims = yLims_raw{i};
    
    if sub0 == 1; ylims = ylims - diff(ylims)/2; end
    
    figure('Position', [100 10 850 9500]); hold on;
    set(gcf, 'renderer', 'Painters')
    
    for sess = sessions
        % Average traces aligned to Motion onset
        sfh1 = subplot(rows,colms,(sess-1)*colms+pSL); hold on
        minY = 0; maxY = 0;
        meanTracesSL_bigSbigT
        xlim(xLimsSL)
        xticks([0 0.2 0.4 0.6]);
        if sess == sessions(end)
            xlabel('Time from motion onset [s]');
        else
            xticklabels({'', '', '', ''})
        end
        if sess == 1
            yPos = ylims(2) + diff(ylims)*0.1 * length(plotSig)/3;
            yPos = ylims(2) + diff(ylims)*0.1 * length(sessions)/6;
            text(-0.02, yPos, '\it Motion onset', 'FontSize', 9)
            
            for l = 1 : 11
                
                if sub0 == 0
                    switch dim
                        case 'MinC'
                            yPos = ylims(2)/2.4-0.03*diff(ylims)*l;
                            xLine = [0.1 0.14];
                        case 'MinI'
                            yPos = ylims(2)/2-0.03*diff(ylims)*l;
                            xLine = [0.1 0.14];
                        otherwise
                            yPos = 1-0.03*l;
                            xLine = [0.03 0.07];
                    end
                else
                    yPos = ylims(2)*0.9-0.03*diff(ylims)*l;
                    xLine = [0.03 0.07];
                end
                line(xLine, yPos*[1 1], 'Color', colorsYeBl_ord{l}, 'LineWidth', 3)
                xPos = xLine(1) + 0.06; legFont = 4;
                if l == 1
                    text(xPos, yPos, 'strong contra', 'FontSize', legFont)
                elseif l == 6
                    text(xPos, yPos, 'weak', 'FontSize', legFont)
                elseif l == 11
                    text(xPos, yPos, 'strong ipsi', 'FontSize', legFont)
                end
            end
            
        end
        ylim(ylims)
        set(gca,'TickDir','out', 'FontSize', fontSize);
        posi = sfh1.Position;
        slWid = sfh1.Position(3); slWid_dim(d) = slWid;
        yposL = ylims(1) + 0.01 * diff(ylims);
        
        if plotGrayLine == 1
            line([0.2 0.5], yposL * [1 1], 'Color', 0.5 * [1 1 1 1], 'LineWidth', 3)
        end
        if colms == 1 || colms == 2
            %ylabel(sigName{i}, 'Fontsize', 9)
            ylabel(['Session ' num2str(sess)], 'Fontsize', 9)
            hYlabel = get(gca, 'YLabel');
            set(hYlabel, 'rotation', 0, 'HorizontalAlignment','right')
        end
        
        % Average traces aligned to Saccade
        if plot_resp_aligned == 1
            sfh2 = subplot(rows,colms,(sess-1)*colms+pSL+1); hold on
            meanTracesRL_bigSbigT
            xlim(xLimsRL)
            xticks([-.2 0]);
            if sess == sessions(end)
                xlabel('Time to saccade [s]');
            else
                xticklabels({'', ''})
            end
            if sess == 1
                yPos = ylims(2) + diff(ylims)*0.1 * length(sessions)/6;
                text(-0.05, yPos, '\it Saccade', 'FontSize', 9)
            end
            
            line([0 0], ylims, 'Color', 'k')
            ylim(ylims)
            set(gca,'TickDir','out', 'FontSize', fontSize);
            rlWid = slWid / abs(diff(xLimsSL)) * abs(diff(xLimsRL));
            posi2 = sfh2.Position;
            posi2(3) = rlWid;
            sfh2.Position = posi2;
            h = gca; h.YAxis.Visible = 'off';
        end
    end
    suptitle(sigName{i})
end
