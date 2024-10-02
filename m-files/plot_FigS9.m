% plot_FigS9


load(fullfile(saveLoc, 'par_sess'))

load(fullfile(saveLoc, 'decoder_excl_SL'))
% Load: 'w_dec', 'pc_dec', 'se_dec', 'excl', 'exclName', 'winCenter','pc_whenC_sign_dec' % Calculated in calc_whatDecoder_fixed_time

load(fullfile(saveLoc, 'decoder_excluded_450'), 'pc450_dec', 'se450_dec') % Calculated in calc_whatDecoder_fixed_time
% Load: 'pc450_dec', 'se450_dec'

rows = 1; colms = 2;
xLimsRL_dec = [-.3 0];
xLimsSL_dec = [0 0.5];

%% Plot Figure 4 a, left: decoder performance aligned to motion onset

figure('Position', [600 400 600 300]); hold on;
plot_excl = [1 12 5 15];
% plot_excl = [9 10 14 11];
clear pltMat pltDat
pltDat = {pc_dec, pc450_dec};
pltColors = {0.5*[1 1 1], [243 158 17]/255, [76 0 153]/255, [48 136 193]/255};
dash = {'-', '--'};
plot_sem = 0;
% -------------------------------
% Stimulus-locked - same time point
sfh1 = subplot(rows, colms, 1); hold on
for plt = 1 : length(pltDat)
    for ee = 1 : length(plot_excl)
        e = plot_excl(ee);
        pltMat = nan(8, size(pltDat{plt}{e, 1}, 2));
        for sess = 1 : 8
            if length(useN_dec{e, sess}) > 0
                pltMat(sess, :) = pltDat{plt}{e, sess};
            end
        end
        
        pltMean = nanmean(pltMat, 1);
        plotCol = pltColors{ee};
        if plot_sem == 0
            plot(winCenter, pltMean, 'LineStyle',dash{plt}, 'Color', plotCol, 'LineWidth', 1)
        else
            
            pltSEM = nanstd(pltMat, [], 1) / sqrt(size(pltMat, 1));
            s = shadedErrorBar(winCenter, pltMean, pltSEM);
            set(s.edge,'LineWidth',0.1,'LineStyle',dash{plt}, 'Color', plotCol)
            s.mainLine.LineWidth = 1.5;
            s.patch.FaceColor = plotCol;
            s.mainLine.Color = plotCol;
        end
        
        
        clear pltMat
    end
end
ylim([0.45 1])
xlim(xLimsSL_dec)
xticks([0 0.2 0.4]);
xlabel('Time from motion onset [s]')
ylabel('Decoding performance')
leg = legend('all neurons', 'only T_{in}', 'no T_{in}', 'no T_{in} or M_{in}', 'Location', 'NorthWest');
set(leg, 'Box', 'off')
set(gca,'TickDir','out', 'FontSize', fontSize);
clear pltDat
grid on
sfh1.Position = [0.1300 0.15 0.35 0.775];
slWid = sfh1.Position(3);


%% Plot Figure 4 a, right: decoder performance aligned to saccade

load(fullfile(saveLoc, 'decoder_excl_RL'),...
    'trl_dec', 'pc_decR', 'se_decR', 'pc450_decR', 'se450_decR', 'excl', 'exclName', 'pc_whenC_sign_decR') % Calculated in calc_whatDecoder_fixed_time

pltDat = {pc_decR, pc450_decR};

sfh2 = subplot(rows, colms, 2); hold on
for plt = 1 : length(pltDat)
    for ee = 1 : length(plot_excl)
        e = plot_excl(ee);
        pltMat = nan(8, size(pltDat{plt}{e, 1}, 2));
        for sess = 1 : 8
            if length(useN_dec{e, sess}) > 0
                pltMat(sess, :) = pltDat{plt}{e, sess};
            end
        end
        
        pltMean = nanmean(pltMat, 1);
        plotCol = pltColors{ee};
        if plot_sem == 0
            plot(trl_dec, pltMean, 'LineStyle',dash{plt}, 'Color', plotCol, 'LineWidth', 1)
        else
            pltMean = nanmean(pltMat, 1);
            pltSEM = nanstd(pltMat, [], 1) / sqrt(size(pltMat, 1));
            s = shadedErrorBar(trl_dec, pltMean, pltSEM);
            set(s.edge,'LineWidth',0.1,'LineStyle','-', 'Color', plotCol)
            s.mainLine.LineWidth = 1.5;
            s.patch.FaceColor = plotCol;
            s.mainLine.Color = plotCol;
        end
        
        clear pltMat
    end
end
xlim(xLimsRL_dec)
ylim([0.45 1])
xticks([-.2 0]);
xlabel('Time to saccade [s]')
grid on
h = gca; h.YAxis.Visible = 'off';
rlWid = slWid / abs(diff(xLimsSL)) * abs(diff(xLimsRL));
sfh2.Position = [0.5703 0.15 rlWid 0.775];
set(gca,'TickDir','out', 'FontSize', fontSize);
clear pltDat







