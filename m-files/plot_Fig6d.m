% plot_Fig6d

% Plot leverage of Min activity on behavior

%% Load data
medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));

%% Plot settings

mediatedDimsA = [3 11 27]; % MinC, MinI, Tin ('for_plots_mediation_norm_to_se')

rows = 1; colms = 2;
sFont = 6;

xLims = [0.2 0.55];
yLimCho = [-3 10];
yLimRT = [-0.4 0.1];

plotSess = [1 : 8];
plot_cols = {[76 0 153]/255, [199 190 0]/255, [41 174 137]/255};
x_shift = 0.005;    % shift of value at mediating time point for visualization purpose

lwR = 1.2; lwM = 0.8;

%% Plot Figure 6 d
figure('Position', [200 10 700 300]); hold on
rr = 0;
for mDim = mediatedDimsA
    
    disp([medA.leverage(mDim).signal])
    
    rr = rr + 1;
    t = medA.leverage(mDim).t;
    plotT = find(t >=0.2 & t <=0.55);
    % compTP = find(round(t, 3) == 0.55);
    plotCol = plot_cols{rr};
    
    % Choice =============================================================
    % Leverage on choice: neuron classes
    subplot(rows, colms, 1); hold on
    
    % Raw
    m = medA.leverage(mDim).choice.unmediated.m; if mDim==60; m = -m; end
    se = medA.leverage(mDim).choice.unmediated.se;
    
    
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.1,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwR;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;

    ylabel('Leverage on choice (slope beta)')
    if rr == length(mediatedDimsA); xlabel('Time from motion onset [s]'); end
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimCho)
    set(gca,'TickDir','out', 'FontSize', sFont);
    
    % RT ==================================================================
    % Correlation with RT
    subplot(rows, colms, 2); hold on
    
    % Raw
    m = medA.leverage(mDim).RT.unmediated.m; if mDim==60; m = -m; end
    se = medA.leverage(mDim).RT.unmediated.se;
    
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.1,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwR;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    ylabel('Correlation with RT (r)')
    if rr == length(mediatedDimsA); xlabel('Time from motion onset [s]'); end
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimRT)
    set(gca,'TickDir','out', 'FontSize', sFont);
    
end


