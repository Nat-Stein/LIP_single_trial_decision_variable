% Plot activity during signal viewing task

i = 0;
for d = 1 : length(plotSig)
    
    dim = plotSig{d};
    i = i + 1;
    
    assignSignal_view
    
     if sub0 == 1; ylims = ylims - diff(ylims)/2; end
     
     % Average traces aligned to Motion onset
     s_view = subplot(rows,colms,(i-1)*colms+pSL); hold on
     minY = 0; maxY = 0;
     meanTracesSLview_bigSbigT
     xlim(xLimView)
     xticks([0 0.2 0.4]);
     %ylabel('Activity [a.u.]')
     if d == length(plotSig)
         xlabel('Time from motion onset [s]');
     else
         xticklabels({'', '', '', ''})
     end
     if d == 1
         yPos = ylims(2) + diff(ylims)*0.1 * length(plotSig)/3;
         text(-0.02, yPos, '\it Motion onset', 'FontSize', 9)
         
         for l = 1 : 11
             
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
    
     if pSL == 1
         ylabel(sigName{i}, 'Fontsize', 9)
         hYlabel = get(gca, 'YLabel');
         set(hYlabel, 'rotation', 0, 'HorizontalAlignment','right')
     end
     
end
