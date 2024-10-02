% Plot bar graphs of mean across session mediation RT only

barCols = {[233 149 53]/255, [76 0 153]/255, [0 153 76]/255, [233 53 149]/255,...
    [0 149 233]/255, [76 0 153]/255, [0 153 76]/255, [233 53 149]/255};
Y = [];
labs = [];
cols = [];
d = 0; 
for dim2 = medDims
    d = d + 1;
    Y = [Y XpercMed_adj{dim, dim2}];
    labs{d} = pltName{dim2};
    if dim == dim2
        cols = [cols [0;0;0]]; 
    else
    cols = [cols barCols{d}']; 
    end
end
cols = cols';
numS = size(Y, 1);
allS = 1 : numS;
stand = allS ~= exS;

if plotBar == 1
    x = [1 1; 2 2];
    b = bar(x, y, 'FaceColor', 'flat', 'FaceAlpha', 0.8);
    for bc = 1 : 2
        b(bc).FaceColor = barColors{bc};
    end
    xEr(:, 1) = x(:, 1)-0.15;
    xEr(:, 2) = x(:, 1)+0.15;
    %er = errorbar(xEr,y,y-ySEM,y+ySEM);
    %er.Color = [0 0 0];
    %er.LineStyle = 'none';
    
elseif plotViolin == 1
    
    violin(Y,'xlabel',labs, 'facealpha', 0.6,...
        'facecolor',cols,'edgecolor','b',...
        'bw',7,...
        'mc',[],...
        'medc','k--') %
    for d2 = 1 : size(Y, 2)
        xPos = d2*ones(1, sum(stand)) + allS(stand) * 0.01 - 0.04;
        plot(xPos, Y(stand, d2), 'o', 'MarkerSize', 10, 'Color', [1 1 1], 'MarkerFaceColor', cols(d2, :))
        plot(d2, Y(exS, d2), 'o', 'MarkerSize', 11, 'Color', [0 0 0], 'MarkerFaceColor', cols(d2, :))
    end
    
elseif plotBox == 1
    
    boxplot(Y,...
        'Width', 0.3, 'BoxStyle', 'outline',... %, 'filled',...
        'ColorGroup', [1 : length(medDims)], 'Colors', [barColors{1}; barColors{2}],...
        'Symbol', '.', 'OutlierSize', 10)
    % 'Labels', {'RT: self','RT: Tin','Choice: self','Choice: Tin'}, 'LabelOrientation', 'horizontal')
end

ax = gca; 
ax.FontSize = fnt; 