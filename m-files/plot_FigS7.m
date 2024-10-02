% plot_FigS7




%% Plot distribution of weights

usePerc = 1 : 2;
colms = 3;

distr_detect = 2;  % Which weights to use for predicted proportion of Tin: 1 = raw weights; 2 = perceptile within session

dimms = [63 20 61 7];


clear ws
figure('Position', [200 10 800 600]); hold on
histbins = [-0.5 : 0.02 : 0.5];
histDims{1} = {[-0.2 : 0.01 : 0.2], histbins, [-0.3 : 0.02 : 0.3], [-0.5 : 0.05 : 0.5]};
histDims_perc = [0 : 0.05 : 1];
histDims{2} = {histDims_perc, histDims_perc, histDims_perc, histDims_perc};

xLims = {[-.2 .2], [-.5 .5], [-.4 .4], [-.5 .5]};
% collies = {[41 174 137]/255, [219 193 22]/255, 0.5 * [1 1 1], [55 189 184]/255};
collies = {[76 0 153]/255, [199 190 0]/255, 0.5 * [1 1 1], [55 189 184]/255};

w_TinC_orig = sigSL.w.TinC;
w_TinI_orig = sigSL.w.TinI;
w_MinC_orig = sigSL.w.MinC;
w_MinI_orig = sigSL.w.MinI;
dd = 0;
for d = 1 : length(dimms)
    dd = dd + 1;
    dim = dimms(d);
    
    w_not = [];
    w_other = [];
    w_TinC = [];
    w_TinI = [];
    w_notTinC = [];
    w_MinC = [];
    w_MinI = [];
    perc_TinC = [];
    perc_TinI = [];
    perc_not = [];
    perc_notTin = [];
    w_allS = [];
    for sess = 1 : 8
        
        w_orig = [];
        if dim == 63
            w_orig = sigSL.w.ramp;
        elseif dim == 20
            w_orig = sigSL.w.PC1;
        elseif dim == 61
            w_orig = sigSL.w.whatD;
        elseif dim == 7
            w_orig = sigSL.w.whenC;
        else
            error('Specify dimension')
        end
        lenW = sum(w_orig{sess}.*w_orig{sess});
        w_orig{sess} = w_orig{sess} / sqrt(lenW);
        
        perc = [];
        for n = 1 : length(w_orig{1, sess})
            a = w_orig{1, sess};
            perc(n) = length(find(a < a(n)))/length(a);
        end
        sz = size(w_orig{1, sess}(w_TinC_orig{1, sess}==0 & w_TinI_orig{1, sess}==0));
        catDim = find(sz ~= 1);
        w_allS = cat(catDim, w_allS, w_orig{1, sess});
        w_not = cat(catDim, w_not, w_orig{1, sess}(w_TinC_orig{1, sess}==0 & w_TinI_orig{1, sess}==0));
        w_other = cat(catDim, w_other, w_orig{1, sess}(w_TinC_orig{1, sess}==0 & w_TinI_orig{1, sess}==0 &...
            w_MinC_orig{1, sess}==0 & w_MinI_orig{1, sess}==0));
        w_notTinC = cat(catDim, w_notTinC, w_orig{1, sess}(w_TinC_orig{1, sess}==0));
        w_TinC = cat(catDim, w_TinC, w_orig{1, sess}(w_TinC_orig{1, sess}~=0));
        w_TinI = cat(catDim, w_TinI, w_orig{1, sess}(w_TinI_orig{1, sess}~=0));
        w_MinC = cat(catDim, w_MinC, w_orig{1, sess}(w_MinC_orig{1, sess}~=0));
        w_MinI = cat(catDim, w_MinI, w_orig{1, sess}(w_MinI_orig{1, sess}~=0));
        
        perc_TinC = cat(2, perc_TinC, perc(w_TinC_orig{1, sess}~=0));
        perc_TinI = cat(2, perc_TinI, perc(w_TinI_orig{1, sess}~=0));
        perc_not = cat(2, perc_not, perc(w_TinC_orig{1, sess}==0  & w_TinI_orig{1, sess}==0));
    end
    
    w_mean(1, dd) = nanmean(w_TinC);
    w_mean(2, dd) = nanmean(w_TinI);
    w_mean(3, dd) = nanmean(w_MinC);
    w_mean(4, dd) = nanmean(w_MinI);
    w_mean(5, dd) = nanmean(w_not);
    w_mean(6, dd) = nanmean(w_other);
    
    
    w_x{dim, 1} = {w_not, w_TinC, w_TinI, w_allS};
    if size(perc_not, 2) < size(perc_not, 1); perc_not=perc_not'; perc_TinC = perc_TinC'; end
    w_x{dim, 2} = {perc_not, perc_TinC, perc_TinI, [perc_not perc_TinC]};
    
    r = 0;
    % ------------------------------------------------------------------
    % Distribution of raw weights
    for m = usePerc
        r = r + 1; subplot(colms, length(dimms), (r-1)*length(dimms) + d); hold on
        thisC = collies{3};
        ddd = w_x{dim, m}{1};
        h1 = histogram(ddd, histDims{m}{d}, 'FaceColor', thisC, 'EdgeColor', thisC);
        for s = 1 : 2
            thisC = collies{s};
            ddd = w_x{dim, m}{s+1};
            h1 = histogram(ddd, histDims{m}{d}, 'FaceColor', thisC, 'EdgeColor', thisC);
        end
        set(gca,'TickDir','out', 'FontSize', sFont);
        xlim([min(histDims{m}{d}) max(histDims{m}{d})]);
        if dim == dimms(1); ylabel('Number of cells'); end
        if m == 1
            xlabel('Weight')
        else
            xlabel('Weight percentile')
        end
    end
    
    % ------------------------------------------------------------------
    % Plot probability of detecting a Tin neuron by applying a
    % threshold on (percentile of) weight
    r = r + 1; subplot(colms, length(dimms), (r-1)*length(dimms) + d); hold on
    for s = 1 : 2
        w_signal = w_x{dim, distr_detect}{s + 1};
        w_ctrl = w_x{dim, distr_detect}{1};
        thisC = collies{s};
        
        if size(w_ctrl, 2) > size(w_ctrl, 1); w_ctrl = w_ctrl'; w_signal = w_signal'; end
        x = [w_ctrl' w_signal'];
        y = [zeros(1, length(w_ctrl)) ones(1, length(w_signal))];
        [mdl]     =  fitglm(x, y,'Distribution', 'binomial','Link','logit');
        x_out = histDims_perc;
        betas = [mdl.Coefficients{1, 1} mdl.Coefficients{2, 1}];
        y_out = glmval(betas',  x_out, 'logit');
        plot(x_out, y_out, 'Color', thisC, 'LineWidth', 1.5);
    end
    ylim([0 1])
    xlim([min(histDims_perc) max(histDims_perc)]);
    set(gca,'TickDir','out', 'FontSize', sFont);
    if dim == dimms(1)
        ylabel('Predicted proportion T_{in}'); 
        leg = legend('T_{in}^{con}', 'T_{in}^{ips}', 'Location', 'NorthWest');
        set(leg, 'Box', 'off')
    end
    xlabel('Weight percentile')
    
end

%% Average weights for TinCon, TinIpsi, MinLeft, MinRight and other

names = {'Ramp', 'PC1', 'When', 'What'};
type = {'TinC', 'TinI', 'MinC', 'MinI', 'no Tin', 'no Tin/Min'};
for tp = 1 : 6
    disp(type{tp})
    for dd = 1 : 4
        disp([names{dd} ': ' num2str(squeeze(w_mean(tp, dd)))])
        
    end
end









