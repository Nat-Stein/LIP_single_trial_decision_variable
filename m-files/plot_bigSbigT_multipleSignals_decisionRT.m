% plot_bigSbigT_multipleSignals_decisionRT


for dim = dimms
    i = i + 1;
     switch dim
        case 1
            sl = sigSL.TinC*1000;
            rl = sigRL.TinC*1000;
            ylims = [0 30]; 
         case 2
             sl = sigSL.TinI;
             rl = sigRL.TinI;
             ylims = [0 20];
         case 7
             sl = sigSL.whenC;
             rl = sigRL.whenC;
             ylims = [0 1.1];
         case 20
             sl = sigSL.PC1;
             rl = sigRL.PC1;
             ylims = [0 1.1];
         case 200
             sl = sigSL.PC1new;
             rl = sigRL.PC1new;
             %ylimsN = [0 1.1];
         case 31
             sl = sigSL.MinC*1000;
             rl = sigRL.MinC*1000;
             ylims = [0 25];
         case 32
             sl = sigSL.MinI*1000;
             rl = sigRL.MinI*1000;
             ylims = [0 30];
         case 61
             sl = sigSL.whatD;
             rl = sigRL.whatD;
             ylims = [-.1 1.1];
         case 63
             sl = sigSL.ramp;
             rl = sigRL.ramp;
             ylims = [-.1 1.1];
         case 70
             sl = sigSL.whenC_tin;
             rl = sigRL.whenC_tin;
             ylims = [0 1.1];
         case 71
             sl = sigSL.whenC_noTin;
             rl = sigRL.whenC_noTin;
             ylims = [0 1.1];
         case 100
             sl = sigSL.MinC*1000 - sigSL.MinI*1000;
             rl = sigRL.MinC*1000 - sigRL.MinI*1000;
             ylims = [-20 20];
         case 101
             sl = cumsum(sigSL.MinC - sigSL.MinI, 2);
             rl = cumsum(sigRL.MinC - sigRL.MinI, 2);
             ylims = [-3 3];
         case 102
             sl = sigSL.MinC - sigSL.MinI;
             rl = sigRL.MinC - sigRL.MinI;
             BL = find(par.tsl >= -0.1 & par.tsl < 0);
             BL_minCI = nanmean(sl(:, BL), 2);
             
             sl = sl - repmat(BL_minCI, 1, size(sl, 2));
             rl = rl - repmat(BL_minCI, 1, size(rl, 2));
             
             sl(:, find(par.tsl<=0.1)) = 0;
             
             sl = cumsum(sl, 2);
             rl = cumsum(rl, 2);
             ylims = [-3 3];
     end
     if sub0 == 1; ylims = ylims - diff(ylims)/2; end
     
     % Stimulus-locked
     sfh1 = subplot(rows,colms,(i-1)*colms+1); hold on
     minY = 0; maxY = 0;
     meanTracesSL_bigSbigT_RT
     xlim(xLimsSL)
     ylim(ylims)
     set(gca,'TickDir','out', 'FontSize', fontSize);
     posi = sfh1.Position;
     slWid = sfh1.Position(3);
     
     yticks(yTicks{i})
     ylabel(sigName{i}, 'Fontsize', 9)
     hYlabel = get(gca, 'YLabel');
     set(hYlabel, 'rotation', 0, 'HorizontalAlignment','right')
     xticks([0 .2 .4 .6]);
     if i == length(plotSig)
         xlabel('Time from motion onset [s]');
     else
         xticklabels({'', '','', ''})
         if i == 1
             yPos = ylims(2) + diff(ylims)*0.1;
             text(-0.04, yPos, '\it Motion onset', 'FontSize', 9)
         end
     end
     
     % Response-locked
     sfh2 = subplot(rows,colms,(i-1)*colms+2); hold on
     meanTracesRL_bigSbigT_RT
     xlim(xLimsRL)
     ylim(ylims)
     set(gca,'TickDir','out', 'FontSize', fontSize);
     rlWid = slWid / abs(diff(xLimsSL)) * abs(diff(xLimsRL));
     posi2 = sfh2.Position;
     posi2(3) = rlWid;
     sfh2.Position = posi2;
     
     xticks([-.2 0]); 
     if dim == dimms(end)
         xlabel('Time to saccade [s]');
     else
         xticklabels({'', ''})
         if dim == dimms(1)
             yPos = ylims(2) + diff(ylims)*0.1;
             text(-0.1, yPos, '\it Saccade', 'FontSize', 9)
         end
     end
     h = gca; h.YAxis.Visible = 'off';
     line([0 0], ylims, 'Color', 'k')
end


