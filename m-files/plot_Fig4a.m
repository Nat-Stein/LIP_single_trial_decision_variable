


load(fullfile(saveLoc, 'par_sess'))

load(fullfile(saveLoc, 'decoder_excl_SL'))
% Load: 'w_dec', 'pc_dec', 'se_dec', 'excl', 'exclName', 'winCenter','pc_whenC_sign_dec' % Calculated in calc_whatDecoder_fixed_time

load(fullfile(saveLoc, 'decoder_excluded_450'), 'pc450_dec', 'se450_dec') % Calculated in calc_whatDecoder_fixed_time
% Load: 'pc450_dec', 'se450_dec'

rows = 1; colms = 2; 
xLimsRL_dec = [-.3 0];
xLimsSL_dec = [0 0.5];

%% Plot Figure 4 a, left: decoder performance aligned to motion onset

figure('Position', [600 400 700 300]); hold on;
e = 1; % Don't exclude any neurons from decoder
clear pltMat pltDat
pltDat = {pc_dec, pc450_dec, pc_whenC_sign_dec};
pltColors = {[243 158 17]/255, [176 28 98]/255, [15 137 103]/255};
% -------------------------------
% Stimulus-locked - same time point
sfh1 = subplot(rows, colms, 1); hold on
for plt = 1 : length(pltDat)
    for sess = 1 : 8
        pltMat(sess, :) = pltDat{plt}{e, sess};
    end
    
    pltMean = nanmean(pltMat, 1);
    pltSEM = nanstd(pltMat, [], 1) / sqrt(size(pltMat, 1));
    s = shadedErrorBar(winCenter, pltMean, pltSEM);
    plotCol = pltColors{plt};
    set(s.edge,'LineWidth',0.1,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = 1.5;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    
    clear pltMat
end
ylim([0.45 1])
xlim(xLimsSL_dec)
xticks([0 0.2 0.4]); 
xlabel('Time from motion onset [s]')
ylabel('Decoding accuracy')
leg = legend('time-dependent', 'fixed training time', 'When-decoder', 'Location', 'NorthWest');
set(leg, 'Box', 'off')
set(gca,'TickDir','out', 'FontSize', fontSize);
clear pltDat
grid on
sfh1.Position = [0.1300 0.15 0.35 0.775];
slWid = sfh1.Position(3);


%% Plot Figure 4 a, right: decoder performance aligned to saccade

load(fullfile(saveLoc, 'decoder_excl_RL'),...
    'trl_dec', 'pc_decR', 'se_decR', 'pc450_decR', 'se450_decR', 'excl', 'exclName', 'pc_whenC_sign_decR') % Calculated in calc_whatDecoder_fixed_time


pltDat = {pc_decR, pc450_decR, pc_whenC_sign_decR}; 

sfh2 = subplot(rows, colms, 2); hold on
for plt = 1 : length(pltDat)
    for sess = 1 : 8
        pltMat(sess, :) = pltDat{plt}{e, sess};
    end
    
    pltMean = nanmean(pltMat, 1);
    pltSEM = nanstd(pltMat, [], 1) / sqrt(size(pltMat, 1));
    s = shadedErrorBar(trl_dec, pltMean, pltSEM);
    plotCol = pltColors{plt};
    set(s.edge,'LineWidth',0.1,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = 1.5;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    clear pltMat
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







