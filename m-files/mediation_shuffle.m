
% run_mediation_all2all


load(fullfile(saveLoc, 'sig_allSessions'),...
    'sigSL')

% Analysis parameters
max_coh = 0.1;
minRT = 0.67;
maxRT = 2;
start_t = 0;
end_t = 0.6;

S = sigSL;
S.MinI_minus_C = S.MinI - S.MinC;
S.Int_MinC_minus_I = cumsum(S.MinC - S.MinI, 2);

session = S.session;
times = S.t;
dt = times(2) - times(1);
choice = S.choice;
RT = S.rt;
coh = S.sig_coh;
cohs = unique(coh)';

numShuffle = 1000;
for i = 1 : numShuffle
    thisShuff = nan(size(session));
    for sess = 1 : 8
        for c = cohs
            trls = find(session == sess & coh == c);
            trls_shuff = trls(randperm(length(trls)));
            thisShuff(trls) = trls_shuff;
        end
    end
    shuff_ord(i, :) = thisShuff;
end
% save(fullfile(saveLoc, 'med_shuff_ord'), 'shuff_ord')

str =  {'MinC',...
    'MinI',...
    'PC1',...
    'TinC',...
    'whatD',...
    'whenC',...
    'ramp',...
    'MinI_minus_C',...
    'Int_MinC_minus_I'};

plotFigs = 0;

% smooth in window
for i = 1 : length(str)
    sm = round(0.05/dt);
    h = ones(sm,1)/sm;
    S.(str{i}) = conv2(1, h, S.(str{i}), 'same');
end

%% time subset
tind = findclose(times, -0.1:0.01:0.75);
for i=1:length(str)
    aux = S.(str{i});
    S.(str{i}) = aux(:,tind);
end
times = times(tind);


%%

flags.norm_to_se = 1;
nsessions =  8;
clear med_shuffle
tic
for sh = 1 : numShuffle
    disp(num2str(sh)); toc
    RT_shuff = RT(shuff_ord(sh, :));
    choice_shuff = choice(shuff_ord(sh, :));
    clear leverage
    
    count = 0;
    for iMethod=1:length(str)
        disp(['Shuffle #' num2str(sh) ', method' num2str(iMethod)]); toc
        sProj = S.(str{iMethod});
        disp(num2str(iMethod));
        for iMediator=1:length(str)
            
            if iMediator~=iMethod
                
                tind = findclose(times, 0.55);
                sProj_mediator = sProj(:,tind);
                sOther_mediator = S.(str{iMediator})(:,tind);
                
                count = count + 1;
                
                %% all sessions together
                
                
                
                %% per session and average
                do_plot = 0;
                clear out out_other
                i = 0;
                for j=1:nsessions
                    I = session==j & abs(coh)<max_coh; % only low coh
                    if ~[all(isnan(to_vec(sProj(I,:)))) || all(isnan(sProj_mediator(I))) || all(isnan(sOther_mediator(I)))]
                        i = i+1;
                        [~,out(i)] = corr_with_RT_choice(choice_shuff(I)==0, RT_shuff(I), coh(I), times, ...
                            sProj(I,:), sProj_mediator(I), minRT,maxRT, start_t, end_t, do_plot,flags);
                        [~,out_other(i)] = corr_with_RT_choice(choice_shuff(I)==0, RT_shuff(I), coh(I), times, ...
                            sProj(I,:), sOther_mediator(I), minRT,maxRT, start_t, end_t,do_plot,flags);
                    end
                end
                
                %% save
                
                leverage(count).signal = str{iMethod};
                leverage(count).mediator = str{iMediator};
                
                leverage(count).t = out(1).tt;
                
                
                x = cat(2, out.rho_RT);
                [leverage(count).RT.unmediated.m , leverage(count).RT.unmediated.se] = averageCorrelation(x);
                
                x = cat(2, out.rho_RT_partial);
                [leverage(count).RT.mediated_self.m, leverage(count).RT.mediated_self.se] = averageCorrelation(x);
                
                x = cat(2, out_other.rho_RT_partial);
                [leverage(count).RT.mediated_other.m, leverage(count).RT.mediated_other.se] = averageCorrelation(x);
                
                leverage(count).choice.unmediated.m = nanmean(cat(2,out.rho_choice),2);
                leverage(count).choice.unmediated.se = stderror(cat(2,out.rho_choice),2);
                leverage(count).choice.mediated_self.m = nanmean(cat(2,out.rho_choice_partial),2);
                leverage(count).choice.mediated_self.se = stderror(cat(2,out.rho_choice_partial),2);
                leverage(count).choice.mediated_other.m = nanmean(cat(2,out_other.rho_choice_partial),2);
                leverage(count).choice.mediated_other.se = stderror(cat(2,out_other.rho_choice_partial),2);
                
                leverage(count).per_session.self_mediated = out;
                leverage(count).per_session.other_mediated = out_other;
                
                med_shuffle.rho_RT(sh, count, :) = leverage(count).RT.unmediated.m;
                med_shuffle.rho_RT_partial(sh, count, :) = leverage(count).RT.mediated_self.m;
                med_shuffle.rho_RT_partial_other(sh, count, :) = leverage(count).RT.mediated_other.m;
                med_shuffle.rho_choice(sh, count, :) = leverage(count).choice.unmediated.m;
                med_shuffle.rho_choice_partial(sh, count, :) = leverage(count).choice.mediated_self.m;
                med_shuffle.rho_choice_partial_other(sh, count, :) = leverage(count).choice.mediated_other.m;
                med_shuffle.t = leverage(count).t;
                med_shuffle.signal{count} = leverage(count).signal;
                med_shuffle.mediator{count} = leverage(count).mediator;
                
            end
        end
        
    end
    
    
end
save(fullfile(saveLoc, 'for_plots_mediation_norm_to_se_shuffle'), 'med_shuffle', 'shuff_ord', 'leverage')



%% Plot settings

mediatedDimsA = [52 19 27]; dimTin = 27;% ramp, PC1, Tin 

rows = length(mediatedDimsA); colms = 2;
sFont = 6;

xLims = [0.2 0.56];
yLimCho = [-6 10];
yLimRT = [-0.4 0.1];

plotSess = [1 : 8];
col_self = [48 136 193]/255;
col_Tin = [237 190 4]/255;
x_shift = 0.005;    % shift of value at mediating time point for visualization purpose

lwR = 1.2; lwM = 0.8;

useShuff = 1 : 500;

%% Plot figure 5
figure('Position', [200 10 700 700]); hold on
rr = 0;
for mDim = mediatedDimsA
    
    %     if mDim ~= dimTin
    %         disp([medA.leverage(mDim).signal ' mediated by self and ' medA.leverage(mDim).mediator])
    %     else
    %         disp([medA.leverage(mDim).signal ' mediated by self'])
    %     end
    
    rr = rr + 1;
    t = med_shuffle.t;
    plotT = find(t >=0.2 & t <=0.5);
    compTP = find(round(t, 3) == 0.55);
    
    % Choice =============================================================
    % Leverage on choice: neuron classes
    subplot(rows, colms, (rr-1)*colms + 1); hold on
    
    % Raw
    
    m = squeeze(nanmean(med_shuffle.rho_choice(useShuff, mDim, :),1)); if mDim==60; m = -m; end
    se = squeeze(nanstd(med_shuffle.rho_choice(useShuff, mDim, :),1));
    %     m = medA.leverage(mDim).choice.unmediated.m; if mDim==60; m = -m; end
    %     se = medA.leverage(mDim).choice.unmediated.se;
    
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
    m = squeeze(nanmean(med_shuffle.rho_choice_partial(useShuff, mDim, :),1)); if mDim==60; m = -m; end
    se = squeeze(nanstd(med_shuffle.rho_choice_partial(useShuff, mDim, :),1));
    %     m = medA.leverage(mDim).choice.mediated_self.m; if mDim==60; m = -m; end
    %     se = medA.leverage(mDim).choice.mediated_self.se;
    
    plotCol = col_self;
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwM;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    % Mediated by Tin
    if mDim ~= dimTin
        m = squeeze(nanmean(med_shuffle.rho_choice_partial_other(useShuff, mDim, :),1)); if mDim==60; m = -m; end
        se = squeeze(nanstd(med_shuffle.rho_choice_partial_other(useShuff, mDim, :),1));
        %         m = medA.leverage(mDim).choice.mediated_other.m; if mDim==60; m = -m; end
        %         se = medA.leverage(mDim).choice.mediated_other.se;
        
        plotCol = col_Tin;
        s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = lwM;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        % Raw Tin late
        %         m = medA.leverage(dimTin).choice.unmediated.m;
        %         se = medA.leverage(dimTin).choice.unmediated.se;
        m = squeeze(nanmean(med_shuffle.rho_choice(useShuff, dimTin, :),1));
        se = squeeze(nanstd(med_shuffle.rho_choice(useShuff, dimTin, :),1));
        
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
    m = squeeze(nanmean(med_shuffle.rho_RT(useShuff, mDim, :),1)); if mDim==60; m = -m; end
    se = squeeze(nanstd(med_shuffle.rho_RT(useShuff, mDim, :),1));
    %     m = medA.leverage(mDim).RT.unmediated.m; if mDim==60; m = -m; end
    %     se = medA.leverage(mDim).RT.unmediated.se;
    
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
    m = squeeze(nanmean(med_shuffle.rho_RT_partial(:, mDim, :),1)); if mDim==60; m = -m; end
    se = squeeze(nanstd(med_shuffle.rho_RT_partial(:, mDim, :),1));
    %     m = medA.leverage(mDim).RT.mediated_self.m; if mDim==60; m = -m; end
    %     se = medA.leverage(mDim).RT.mediated_self.se;
    
    plotCol = col_self;
    s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
    set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
    s.mainLine.LineWidth = lwM;
    s.patch.FaceColor = plotCol;
    s.mainLine.Color = plotCol;
    
    % Mediated by Tin
    if mDim ~= dimTin
        m = squeeze(nanmean(med_shuffle.rho_RT_partial_other(useShuff, mDim, :),1)); if mDim==60; m = -m; end
        se = squeeze(nanstd(med_shuffle.rho_RT_partial_other(useShuff, mDim, :),1));
        %         m = medA.leverage(mDim).RT.mediated_other.m; if mDim==60; m = -m; end
        %         se = medA.leverage(mDim).RT.mediated_other.se;
        
        plotCol = col_Tin;
        s = shadedErrorBar(t(plotT), m(plotT), se(plotT));
        set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
        s.mainLine.LineWidth = lwM;
        s.patch.FaceColor = plotCol;
        s.mainLine.Color = plotCol;
        
        % Raw Tin late
        %         m = medA.leverage(dimTin).RT.unmediated.m;
        %         se = medA.leverage(dimTin).RT.unmediated.se;
        m = squeeze(nanmean(med_shuffle.rho_RT(useShuff, dimTin, :),1));
        se = squeeze(nanstd(med_shuffle.rho_RT(useShuff, dimTin, :),1));
        
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




















