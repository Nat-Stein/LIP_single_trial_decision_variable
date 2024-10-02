%% Plot SL and RL coherence for all sessions combined

dimms = [1 2 5 31 32 55 56]; rows = length(dimms); % 33
sc = 1; 

plotFilt = ones(1, 50)/50;
xLimsSL = [0 0.5];
xLimsRL = [-0.5 0];
figure('Position', [600 100 800 900]); hold on; rows = 8; i = 0;
for sess = 1 : 8
    
    
    for dim = dimms
        i = i + 1;  
        subplot(rows, 2, (i-1)*2+1); hold on
        minY = 0; maxY = 0;
        % meanTracesSL_Fig4
        meanTracesSL_bigSbigT
        xlim(xLimsSL)
        ylabel('Activity [a.u.]')
        set(gca,'TickDir','out', 'FontSize', 10);
        title(dimName{dim})
        
        subplot(rows, 2, (i-1)*2+2); hold on
        meanTracesRL_plot
        xlim(xLimsRL)
        ylabel('Activity [a.u.]')
        set(gca,'TickDir','out', 'FontSize', 10);
        
        if minY == maxY; maxY = maxY + 0.001; end
        
        for plt = 1 : 2
            subplot(rows, 2, (i-1)*2+plt); hold on
            %ylim([minY ceil(maxY/5)*5])
            ylim(1.05 * [minY maxY])
        end
        
    end
    %suptitle(par_sess{sess}.date)
    
end