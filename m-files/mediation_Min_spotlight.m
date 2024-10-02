% plot_Figure5

% Plot behavioral correlations based on Ariel's calculations

%% Load data
% medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));
% medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se_test'));
medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se_MinCI_BL'));

% newest PC1 (Dec 20)

%% Plot settings

% mediatedDimsA = [3 11 60 68]; % MinC, MinI, MinI-MinC, integrated MinC-MinI ('for_plots_mediation_norm_to_se_test')
% titles = {'M^{con}_{in}', 'M^{ips}_{in}', 'M^{con}_{in}-M^{ips}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}'};
% dimTin = 24;

xLims = [0.2 0.56]; 
yLimCho = [-3 10];
yLimRT = [-0.4 0.15];

plotSess = [1 : 8];
col_self = [48 136 193]/255;
col_Tin = [237 190 4]/255;
x_shift = 0.005; % shift of value at mediating time point for visualization purpose
sFont = 6;

lwR = 1.2; lwM = 0.8; 


%% Plot only raw values
rows = 1;

mediatedDimsA = [25 60 3 68 11]; % TinC, MinC, MinI, MinI-MinC, integrated MinC-MinI ('for_plots_mediation_norm_to_se_test')
titles = {'T^{con}_{in}', 'M^{con}_{in}', 'M^{ips}_{in}', 'M^{con}_{in}-M^{ips}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}'};

mediatedDimsA = [25 60 3 68 11]; % TinC, MinC, MinI, MinI-MinC, integrated MinC-MinI ('for_plots_mediation_norm_to_se_test')
titles = {'T^{con}_{in}', 'M^{con}_{in}-M^{ips}_{in}', 'M^{con}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}', 'M^{ips}_{in}'};

mediatedDimsA = [3 12 67 76 85]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'M^{con}_{in}', 'M^{ips}_{in}', 'M^{con}_{in}-M^{ips}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}, BL-corr'};
dimTin = 30;

mediatedDimsA = [76 85 91]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'PC1 Min'};
dimTin = 30;

plotCols = {[199 190 0]/255, [93 140 232]/255, [41 174 137]/255, [174 41 98]/255, [76 0 153]/255};

% After addind minReg
mediatedDimsA = [95 105 115]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'Min regression', 'integ Min regression'};

mediatedDimsA = [80 105 115]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'M^{con}_{in}-M^{ips}_{in}', 'Min regression', 'integ Min regression'};
dimTin = 30;

mediatedDimsA = [100 115]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'integ M^{con}_{in}-M^{ips}_{in}, BL-corr, norm'};
dimTin = 30;

mediatedDimsA = [100 40]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'T^{con}_{in}'};
dimTin = 30;


mediatedDimsA = [100 10 40]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'M^{con}_{in}', 'T^{con}_{in}'};
dimTin = 30;


 yLimCho = [-3 10]; yLimRT = [-.4 .15];
%yLimCho = [-3 6]; yLimRT = [-.15 .15];

% plotCols = {[41 174 137]/255, [93 140 232]/255, [76 0 153]/255, [174 41 98]/255, [199 190 0]/255};
legs = []; leg_hand = [];
figure('Position', [200 10 700 250]); hold on
rr = 0;
for mDim = mediatedDimsA
    
    disp([medA.leverage(mDim).signal ' mediated by ' medA.leverage(mDim).mediator])
    
    rr = rr + 1;
    t = medA.leverage(mDim).t;
    plotT = find(t >=-0.2 & t <=0.5);
    compTP = find(round(t, 3) == 0.55);
    
    % Choice =============================================================
    % Leverage on choice: neuron classes
    subplot(rows, colms, 1); hold on
    
    % Raw
    m = medA.leverage(mDim).choice.unmediated.m; % if mDim==60; m = -m; end
    se = medA.leverage(mDim).choice.unmediated.se;
    
    plotCol = plotCols{rr};
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwR;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    plot(t(plotT), m(plotT), 'Color', plotCol, 'LineWidth', 2);
    
    
    pltMeanRaw = squeeze(m(compTP));
    pltSERaw = squeeze(se(compTP));
    if mDim ~= dimTin
        xpos = t(compTP) + x_shift;
    else
        xpos = t(compTP) - x_shift;
    end
    plot(xpos, pltMeanRaw, 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
    line(xpos * [1 1], pltMeanRaw + [-1 1]*pltSERaw, 'Color', plotCol)
 
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
    
    plotCol = plotCols{rr};
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
    
    leg_hand(rr) = plot(t(plotT), m(plotT), 'Color', plotCol, 'LineWidth', 2);
    legs{rr} = titles{rr};
    
    ylabel('Correlation with RT (r)')
    if rr == length(mediatedDimsA); xlabel('Time from motion onset [s]'); end
    line(xLims, [0 0], 'Color', 'k')
    xlim(xLims)
    ylim(yLimRT)
    set(gca,'TickDir','out', 'FontSize', sFont);
    
    tt = medA.leverage.t;
    tpts = find(tt >= 0.2 & tt <= 0.55);
    x = tt(tpts);
    y = medA.leverage(mDim).choice.unmediated.m(tpts)';
    mdl = fitlm(x,y);
    co = table2array(mdl.Coefficients);
    slope = co(2, 1); tval = co(2, 3); pval = co(2, 4);
    disp([medA.leverage(mDim).signal ' on choice: '])
    disp(['Slope = ' num2str(slope)])
    disp(['t = ' num2str(tval) ', p = ' num2str(pval)])
    % Weighted regression
    

    
end
legend(leg_hand, legs)


mediatedDimsA = [100 10 20 38]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'M^{con}_{in}', 'T^{con}_{in}'};
for mDim = mediatedDimsA

    tt = medA.leverage.t;
    tpts = find(tt >= 0.2 & tt <= 0.55);
    clear slope_coeff_sess
    for sess = 1 : length(medA.leverage(mDim).per_session.self_mediated)
        
        ch_rho = medA.leverage(mDim).per_session.self_mediated(sess).rho_choice;
        ch_rho_se = medA.leverage(mDim).per_session.self_mediated(sess).norm_choice;
        
        w = 1./ ch_rho_se(tpts).^2;
        
        x = tt(tpts)';
        y = ch_rho(tpts);
        
        tbl = table(x, y, w, 'VariableNames', {'x', 'y', 'Weights'});
        lm = fitlm(tbl, 'y~x', 'Weights', 'Weights');
        
        all_lm{sess} = lm;
        lm_out = table2array(lm.Coefficients);
        slope_coeff_sess(sess, :, :) = lm_out;
    end
    
    clear all_slopes all_p
    all_slopes = squeeze(slope_coeff_sess(:, 2, 1));
    all_p = squeeze(slope_coeff_sess(:, 2, 4));
    disp('===============')
    disp(medA.leverage(mDim).signal)
    disp(['Highest p = ' num2str(max(all_p))])
    [u, i] = ttest(all_slopes);
    disp(['Mean slopes p = ' num2str(mean(all_slopes))])
    disp(['T-test on slopes p = ' num2str(i)])
end





%% Plot figure 5 equivalent

mediatedDimsA = [3 12 67 76 85]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'M^{con}_{in}', 'M^{ips}_{in}', 'M^{con}_{in}-M^{ips}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}', 'integ M^{con}_{in}-M^{ips}_{in}, BL-corr'};
dimTin = 30;


mediatedDimsA = [103 3 40]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'M^{con}_{in}', 'T^{con}_{in}'};
dimTin = 40;

mediatedDimsA = [103 25]; % 'for_plots_mediation_norm_to_se_MinCI_BL'
titles = {'integ M^{con}_{in}-M^{ips}_{in}, BL-corr', 'PC1'};
dimTin = 40;
rows = length(mediatedDimsA); colms = 2;



figure('Position', [200 10 700 700]); hold on
rr = 0;
for mDim = mediatedDimsA
    
    disp([medA.leverage(mDim).signal ' mediated by ' medA.leverage(mDim).mediator])
    
    rr = rr + 1;
    t = medA.leverage(mDim).t;
    plotT = find(t >=0.2 & t <=0.5);
    compTP = find(round(t, 3) == 0.55);
    
    % Choice =============================================================
    % Leverage on choice: neuron classes
    subplot(rows, colms, (rr-1)*colms + 1); hold on
    title(titles{rr})
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
    title(titles{rr})
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


