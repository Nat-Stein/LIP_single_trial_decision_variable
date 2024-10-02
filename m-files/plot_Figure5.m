% plot_Figure5

% Plot leverage of Sramp, Spc1 and STin on behavior

%% Load data
medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));

%% Plot settings

mediatedDimsA = [52 19 27]; dimTin = 27;% ramp, PC1, Tin ('for_plots_mediation_norm_to_se')
% mediatedDimsA = [3 11 60 68]; % MinC, MinI, MinI-MinC, integrated MinC-MinI ('for_plots_mediation_norm_to_se')

rows = length(mediatedDimsA); colms = 2;
sFont = 6;

xLims = [0.2 0.56]; 
yLimCho = [-2 10];
yLimRT = [-0.4 0.1];

plotSess = [1 : 8];
col_self = [48 136 193]/255;
col_Tin = [237 190 4]/255;
x_shift = 0.005;    % shift of value at mediating time point for visualization purpose

lwR = 1.2; lwM = 0.8; 

%% Plot figure 5
figure('Position', [200 10 700 700]); hold on
rr = 0;
for mDim = mediatedDimsA
    
    if mDim ~= dimTin
        disp([medA.leverage(mDim).signal ' mediated by self and ' medA.leverage(mDim).mediator])
    else
        disp([medA.leverage(mDim).signal ' mediated by self'])
    end
    
    rr = rr + 1;
    t = medA.leverage(mDim).t;
    plotT = find(t >=0.2 & t <=0.5);
    compTP = find(round(t, 3) == 0.55);
    
    % Choice =============================================================
    % Leverage on choice: neuron classes
    subplot(rows, colms, (rr-1)*colms + 1); hold on
    
    % Raw
    m = medA.leverage(mDim).choice.unmediated.m; if mDim==60; m = -m; end
    se = medA.leverage(mDim).choice.unmediated.se;
    
    plotCol = 'k';
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwR;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    pltMeanRaw = squeeze(m(compTP));
    pltSERaw = squeeze(se(compTP));
    if mDim ~= dimTin
        xpos = t(compTP) + x_shift;
    else
        xpos = t(compTP) - x_shift;
    end
    plot(xpos, pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
    line(xpos * [1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    
    % Mediated by self
    m = medA.leverage(mDim).choice.mediated_self.m; if mDim==60; m = -m; end
    se = medA.leverage(mDim).choice.mediated_self.se;
    
    plotCol = col_self;
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwM;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    % Mediated by Tin
    if mDim ~= dimTin
        m = medA.leverage(mDim).choice.mediated_other.m; if mDim==60; m = -m; end
        se = medA.leverage(mDim).choice.mediated_other.se;
        
        plotCol = col_Tin;
        s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = lwM;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        % Raw Tin late
        m = medA.leverage(dimTin).choice.unmediated.m;
        se = medA.leverage(dimTin).choice.unmediated.se;
        
        pltMeanRaw = squeeze(m(compTP));
        pltSERaw = squeeze(se(compTP));
        xpos = t(compTP) - x_shift;
        line(xpos*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'LineStyle', '--', 'Color', [0 0 0])
        plot(xpos, pltMeanRaw, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [1 1 1])
    end
    ylabel('Leverage on choice (slope beta)')
    if rr == length(mediatedDimsA); xlabel('Time from motion onset [s]'); end
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimCho)
    set(gca,'TickDir','out', 'FontSize', sFont);
    
    % RT ==================================================================
    % Correlation with RT
    subplot(rows, colms, (rr-1)*colms + 2); hold on
    
    % Raw
    m = medA.leverage(mDim).RT.unmediated.m; if mDim==60; m = -m; end
    se = medA.leverage(mDim).RT.unmediated.se;
    
    plotCol = 'k';
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwR;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    pltMeanRaw = squeeze(m(compTP));
    pltSERaw = squeeze(se(compTP));
    if mDim ~= dimTin
        xpos = t(compTP) + x_shift;
    else
        xpos = t(compTP) - x_shift;
    end
    plot(xpos, pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
    line(xpos * [1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
    
    % Mediated by self
    m = medA.leverage(mDim).RT.mediated_self.m; if mDim==60; m = -m; end
    se = medA.leverage(mDim).RT.mediated_self.se;
    
    plotCol = col_self;
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwM;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    % Mediated by Tin
    if mDim ~= dimTin
        m = medA.leverage(mDim).RT.mediated_other.m; if mDim==60; m = -m; end
        se = medA.leverage(mDim).RT.mediated_other.se;
        
        plotCol = col_Tin;
        s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = lwM;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        % Raw Tin late
        m = medA.leverage(dimTin).RT.unmediated.m;
        se = medA.leverage(dimTin).RT.unmediated.se;
        
        pltMeanRaw = squeeze(m(compTP));
        pltSERaw = squeeze(se(compTP));
        xpos = t(compTP) - x_shift;
        line(xpos*[1 1], pltMeanRaw + [-1 1]*pltSERaw, 'LineStyle', '--', 'Color', [0 0 0])
        plot(xpos, pltMeanRaw, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [1 1 1])
    end
    
    ylabel('Correlation with RT (r)')
    if rr == length(mediatedDimsA); xlabel('Time from motion onset [s]'); end
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimRT)
    set(gca,'TickDir','out', 'FontSize', sFont);
    
end


