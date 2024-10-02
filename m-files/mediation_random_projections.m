t = round(medA.leverage(1).t, 3);
plotT = find(t >=0.2 & t <=0.5);
compTP = find(round(t, 3) == 0.55);

proj_type = 'PC1weights'; % 'PC1weights', 'uniform', 'randn', 'gauss'

for sess = 1 : 8
    
    clear spk par dat_raw
    load(fullfile(saveLoc, ['spk_S' num2str(sess)]), 'spk', 'par')
    
    par.tsl = round(par.tsl, 3);
    dt = round(par.tsl(2)-par.tsl(1), 3);
    sm = round(0.05/dt);
    h = ones(sm,1)/sm;
    useT = find(par.tsl >= start_t & par.tsl <= end_t);
    t_med = round(par.tsl(useT), 3);
    
    % behavioral indicators for all trials
    choice = par.beh.cho_trg';
    RT = par.beh.rt';
    coh = par.beh.sig_coh';
    
    
    % Subsample to extract only relevant time points
    dat_raw = spk.SLnan(:, :, useT);
    
    load(fullfile(saveLoc, ['PC1_w_S' num2str(sess)]))
    
    for i = 1 : 1000
        
        if rem(i, 50) == 0; disp(['Session ' num2str(sess) ', rep ' num2str(i)]); end
        
        clear dat_proj dat_conv
        
        switch proj_type
            case 'PC1weights'
                % Generate random directions by shuffling weights of PC1
                idx = randperm(length(w));
                w_rand = w(idx);
            case 'uniform'
                % Generste random directions by assigning weights dran from a
                % uniform distribution
                w_rand = rand(1, length(w))-0.5;
            case 'gauss'
                w_rand = randn(1, length(w));
            otherwise
                error('Method does not exist')
        end
        
        w_rand = w_rand / sum(abs(w_rand));
        
        randProj{sess}(i, :) = w_rand;
        
        % Project data onto random vector
        dat_proj = squeeze(mmx('mult', w_rand, dat_raw));
        
        dat_conv = conv2(1, h, dat_proj, 'same');
        
        % Select trials and time points
        tpts = find(ismember(t_med, t));
        I = find(abs(coh)<max_coh); % only low coh
        
        sProj = dat_conv(:, tpts);
        
        sProj_mediator = dat_conv(:, t_med == 0.55);
        
        clear out
        [~,out] = corr_with_RT_choice(choice(I)==0, RT(I), coh(I), t_med(tpts), ...
            sProj(I,:), sProj_mediator(I), minRT, maxRT, start_t, end_t, do_plot,flags);
        
        % Am I taking the right value here???
        med_randProj.rho_RT(sess, i, :) = out.rho_RT;
        med_randProj.rho_RT_partial(sess, i, :) = out.rho_RT_partial;
        med_randProj.rho_choice(sess, i, :) = out.rho_choice;
        med_randProj.rho_choice_partial(sess, i, :) = out.rho_choice_partial;
        
    end
    
end

save(fullfile(saveLoc, ['mediation_randProj_' proj_type]), 'med_randProj', 't_med', 'tpts', 'randProj')


tt = t_med(tpts);
tt_med = find(tt==0.55);
% Plot
figure('Position', [200 10 700 250]); hold on
% Choice =============================================================
% Leverage on choice: neuron classes
subplot(1, 2, 1); hold on
% Raw choice
raw_choice_mean = squeeze(nanmean(med_randProj.rho_choice, 1));
m = nanmean(raw_choice_mean, 1);
%se = nanstd(raw_choice_mean, 1)/sqrt(size(raw_choice_mean, 1));
se = nanstd(raw_choice_mean, 1);
plotCol = 'k';
s = shadedErrorBar(tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(tt(tt_med), m(tt_med), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(tt(tt_med)*[1 1], m(tt_med) + se(tt_med) * [-1 1], 'Color', plotCol)
% Mediated by self
part_choice_mean = squeeze(nanmean(med_randProj.rho_choice_partial, 1));
m = nanmean(part_choice_mean, 1);
% se = nanstd(part_choice_mean, 1)/sqrt(size(part_choice_mean, 1));
se = nanstd(part_choice_mean, 1);
plotCol = col_self;
s = shadedErrorBar(tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwM;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
% Formatting
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim([-2 10])
xlabel('Time from motion onset [s]');
ylabel('Leverage on choice (slope beta)')
set(gca,'TickDir','out', 'FontSize', sFont);
% RT ==================================================================
% Correlation with RT
subplot(1, 2, 2); hold on
plotCol = [0 0 0];
% Raw choice
raw_RT_mean = squeeze(nanmean(med_randProj.rho_RT, 1));
m = nanmean(raw_RT_mean, 1);
%se = nanstd(raw_RT_mean, 1)/sqrt(8);
se = nanstd(raw_RT_mean, 1);
plotCol = 'k';
s = shadedErrorBar(tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(tt(tt_med), m(tt_med), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(tt(tt_med)*[1 1], m(tt_med) + se(tt_med) * [-1 1], 'Color', plotCol)
% Mediated by self
part_RT_mean = squeeze(nanmean(med_randProj.rho_RT_partial, 1));
m = nanmean(part_RT_mean, 1);
% se = nanstd(part_RT_mean, 1)/sqrt(size(part_RT_mean, 1));
se = nanstd(part_RT_mean, 1);
plotCol = col_self;
s = shadedErrorBar(tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwM;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
% Formatting
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim([-0.4 0.1])
xlabel('Time from motion onset [s]');
ylabel('Correlation with RT (r)')
set(gca,'TickDir','out', 'FontSize', sFont);




