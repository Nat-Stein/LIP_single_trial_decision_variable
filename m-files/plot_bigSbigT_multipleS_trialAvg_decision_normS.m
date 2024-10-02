% newer version is called: plot_bigSbigT_multipleS_trialAvg_decision

pSL = colms - 1;

for dim = dimms
    i = i + 1;
    clear sl rl ylims
     switch dim
        case 1
            sl = sigSL.TinC*1000;
            rl = sigRL.TinC*1000;
            ylims = [0 30]; 
            if sub0 == 1
                ylims = [0 20]; 
            end
         case 2
             sl = sigSL.TinI;
             rl = sigRL.TinI;
             ylims = [0 20];
         case 7
             sl = normS.sigSL.whenC;
             rl = normS.sigRL.whenC;
             ylims = [0 1.1];
         case 20
             sl = normS.sigSL.PC1;
             rl = normS.sigRL.PC1;
             ylims = [0 1.1];
             if sub0 == 1
                ylims = [0 0.6]; 
             end
         case 200
             sl = sigSL.PC1new;
             rl = sigRL.PC1new;
             ylims = [-0.5 1];
         case 31
             sl = sigSL.MinC;
             rl = sigRL.MinC;
             ylims = [0 25];
         case 32
             sl = sigSL.MinI;
             rl = sigRL.MinI;
             ylims = [0 30];
         case 61
             sl = normS.sigSL.whatD;
             rl = normS.sigRL.whatD;
             ylims = [-.1 1.1];
         case 63
             sl = normS.sigSL.ramp;
             rl = normS.sigRL.ramp;
             ylims = [-.1 1.1];
             if sub0 == 1
                ylims = [0 0.6]; 
             end
         case 70
             sl = normS.sigSL.whenC_tin; 
             rl = normS.sigRL.whenC_tin;
             ylims = [0 1.1];
         case 71
             sl = normS.sigSL.whenC_noTin; 
             rl = normS.sigRL.whenC_noTin;
             ylims = [0 1.1];
     end
     if sub0 == 1; ylims = ylims - diff(ylims)/2; end
     
     % Average traces aligned to Motion onset
     sfh1 = subplot(rows,colms,(i-1)*colms+pSL); hold on
     minY = 0; maxY = 0;
     meanTracesSL_bigSbigT
     xlim(xLimsSL)
     xticks([0 0.2 0.4 0.6]);
     %ylabel('Activity [a.u.]')
     if dim == dimms(end)
         xlabel('Time from motion onset [s]');
     else
         xticklabels({'', '', '', ''})
     end
     if dim == dimms(1)
         yPos = ylims(2) + diff(ylims)*0.1 * length(dimms)/3;
         text(-0.02, yPos, '\it Motion onset', 'FontSize', 9)
         % hLg = legend({'strong contra', '', '', '', '', 'weak', '', '', '', '', 'strong ipsi'}, 'Location', 'NorthWest', 'FontSize', 4);
         
         for l = 1 : 11
             yPos = 1-0.03*l;
             line([0.03 0.07], yPos*[1 1], 'Color', colorsYeBl_ord{l}, 'LineWidth', 3)
             xPos = 0.09; legFont = 4;
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
     slWid = sfh1.Position(3);
     yposL = ylims(1) + 0.01 * diff(ylims);
     line([0.2 0.5], yposL * [1 1], 'Color', 0.5 * [1 1 1 1], 'LineWidth', 3)
     if colms == 2
         ylabel(sigName{i}, 'Fontsize', 9)
         hYlabel = get(gca, 'YLabel');
         set(hYlabel, 'rotation', 0, 'HorizontalAlignment','right')
     end
     
     % Average traces aligned to Saccade
     sfh2 = subplot(rows,colms,(i-1)*colms+pSL+1); hold on
     meanTracesRL_bigSbigT
     xlim(xLimsRL)
     xticks([-.2 0]);
     if dim == dimms(end)
         xlabel('Time to saccade [s]');
     else
         xticklabels({'', ''})
     end
     if dim == dimms(1)
         yPos = ylims(2) + diff(ylims)*0.1 * length(dimms)/3;
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



