
load(fullfile(saveLoc, 'sig_allSessions'))
load(fullfile(saveLoc, 'par_sess'))

dimms = [31 32 100 101 102]; 
plotSig = {'MinC', 'MinI'};
sigName = {'S^{left}_{Min} [sp/s]', 'S^{right}_{Min} [sp/s]', 'S^{left}_{Min} - S^{right}_{Min} [sp/s]',...
    'Int(S^{left}_{Min} - S^{right}_{Min})',...
    'Int(S^{left}_{Min} - S^{right}_{Min}, BLcorr)'};

par = par_sess{1};
yTicks = {[0 10 20 30], [0 10 20 30], [-20 -10 0 10 20], [-2 0 2], [-2 0 2]}; 
% cutoff = 0.7; ls = '-'; lw = 1;
% smooth_win_t = 0.025;
rows = length(dimms); colms = 2; 

saveFig = 0;
xLimsSL = [-.2 .6];

sub0 = 0; plotGrayLine = 0; 
plot_resp_aligned = 1;

clear qRT
nQuant = 3; 
qRT = quantile(sigSL.rt, nQuant); qRT = [0 qRT max(sigSL.rt)];
rt_quant = nan(size(sigSL.rt)); 
for qu = 1 : length(qRT)-1
    qq = find((sigSL.rt >= qRT(qu) & sigSL.rt <= qRT(qu+1)) & sigSL.choice == 1);
    rt_quant(qq) = qu; 
    qq = find((sigSL.rt >= qRT(qu) & sigSL.rt <= qRT(qu+1)) & sigSL.choice == 0);
    rt_quant(qq) = qu+nQuant+1; 
end

% nQuant = 3; 
% trials = find(sigSL.sig_coh == 0);
% qRT0 = quantile(sigSL.rt(trials), nQuant); qRT0 = [0 qRT0 max(sigSL.rt(trials))];
% rt_quant = nan(size(sigSL.rt)); 
% for qu = 1 : length(qRT0)-1
%     qq = find((sigSL.rt >= qRT0(qu) & sigSL.rt <= qRT0(qu+1)) & sigSL.choice == 1);
%     rt_quant(qq) = qu; 
%     qq = find((sigSL.rt >= qRT0(qu) & sigSL.rt <= qRT0(qu+1)) & sigSL.choice == 0);
%     rt_quant(qq) = qu+nQuant+1; 
% end


% Plot motion discrimination task split by RT percentile
figure('Position', [200 10 600 600]); hold on;  i = 0;
set(gcf, 'renderer', 'Painters')
plot_bigSbigT_multipleSignals_decisionRT


%% Min control analyses

trials = find(sigSL.sig_coh == -0.032);
rt = sigSL.rt(trials);
cho = sigSL.choice(trials);
medRT = median(rt);

trSlow = trials(find(rt > medRT));
trFast = trials(find(rt <= medRT));

TRLs = {find(rt <= medRT & cho == 0),...
    find(rt > medRT & cho == 0);...
    find(rt <= medRT & cho == 1),...
    find(rt > medRT & cho == 1)};
RT_cols =  {[117 13 85]/255, [243 162 193]/255;...
    [57 98 18]/255, [190 217 56]/255};

% colors =  [[29 42 16]/255; [57 98 18]/255; [106 186 0]/255; [190 217 56]/255; [203 231 160]/255;...
%     [104 3 50]/255; [117 13 85]/255; [213 65 123]/255; [243 162 193]/255; [239 207 214]/255];

rows = 4; cols = 2; 
figure; hold on
for c = 1 : 2
    for sf = 1 : 2
        trls = trials(TRLs{sf, c});
        subplot(rows, cols, 1); hold on; title('MinC')
        plt = conv(nanmean(sigSL.MinC(trls, :)), ones(1, 25), 'same') * 100;
        plt(par.tsl > median(sigSL.rt(trls))) = nan;
        plot(par.tsl, plt, 'Color', RT_cols{sf, c})
        subplot(rows, cols, 3); hold on; title('MinC, BL subtracted')
        plt = plt - nanmean(plt(par.tsl >= -0.1 & par.tsl < 0));
        plot(par.tsl, plt, 'Color', RT_cols{sf, c})
        
        subplot(rows, cols, 2); hold on; title('MinI')
        plt = conv(nanmean(sigSL.MinI(trls, :)), ones(1, 25), 'same') * 100;
        plt(par.tsl > median(sigSL.rt(trls))) = nan;
        plot(par.tsl, plt, 'Color', RT_cols{sf, c})
        subplot(rows, cols, 4); hold on; title('MinI, BL subtracted')
        plt = plt - nanmean(plt(par.tsl >= -0.1 & par.tsl < 0));
        plot(par.tsl, plt, 'Color', RT_cols{sf, c})
        
        subplot(rows, cols, 5); hold on; title('MinC-MinI')
        plt = conv(nanmean(sigSL.MinC(trls, :)-sigSL.MinI(trls, :)), ones(1, 25), 'same') * 100;
        plt(par.tsl > median(sigSL.rt(trls))) = nan;
        plot(par.tsl, plt, 'Color', RT_cols{sf, c})
       
        subplot(rows, cols, 7); hold on; title('Int(MinC-MinI)')
        plt1 = plt; plt1(par.tsl<=0) = 0;
        pltInt = cumsum(plt1);
        plot(par.tsl, pltInt, 'Color', RT_cols{sf, c})
        
        subplot(rows, cols, 6); hold on; title('MinC-MinI, BL subtr')
        plt = plt - nanmean(plt(par.tsl >= -0.1 & par.tsl < 0));
        plot(par.tsl, plt, 'Color', RT_cols{sf, c})
        
        subplot(rows, cols, 8); hold on; title('Int(MinC-MinI), BL subtr')
        plt1 = plt; plt1(par.tsl<=0) = 0;
        pltInt = cumsum(plt1);
        plot(par.tsl, pltInt, 'Color', RT_cols{sf, c})
        
    end
end
for p = 1 : 8
    subplot(rows, cols, p); hold on; 
    ylabel('Activity [sp/s]')
    xlabel('Time from motion onset [s]')
    xlim([-0.2 1.5])
end

    
    % Slow trials
    trls = trials(find(rt > medRT & cho == c));
    subplot(1, 2, 1); hold on
    plt = conv(nanmean(sigSL.MinC(trls, :)), ones(1, 25), 'same');
    plt(par.tsl > median(sigSL.rt(trls))) = nan;
    plot(par.tsl, plt)
    
    subplot(1, 2, 2); hold on
    plt = conv(nanmean(sigSL.MinI(trls, :)), ones(1, 25), 'same');
    plt(par.tsl > median(sigSL.rt(trls))) = nan;
    plot(par.tsl, plt)
    
    % Fast trials
    trls = trials(find(rt <= medRT & cho == c));
    subplot(1, 2, 1); hold on
    plt = conv(nanmean(sigSL.MinC(trls, :)), ones(1, 25), 'same');
    plt(par.tsl > median(sigSL.rt(trls))) = nan;
    plot(par.tsl, plt)
    
    subplot(1, 2, 2); hold on
    plt = conv(nanmean(sigSL.MinI(trls, :)), ones(1, 25), 'same');
    plt(par.tsl > median(sigSL.rt(trls))) = nan;
    plot(par.tsl, plt)
    
end
legend('slow left', 'fast left', 'slow right', 'fast right')
















