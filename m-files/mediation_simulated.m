

% Set up format of simulation data to match mediation analysis

medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));
sim_dat = load(fullfile(saveLoc, 'race_model_simulation'));
load(fullfile(saveLoc, 'par_sess'))

% Add time point 0.2s (0)
sim_dat.diffusionPaths = [zeros(60000, 1) sim_dat.diffusionPaths];
sim_dat.tSim = [0 sim_dat.tSim];

start_t = sim_dat.tSim(1) + 0.2;
end_t = 0.6;

% Parameters of mediation analysis
max_coh = 0.1;
minRT = 0.67;
maxRT = 2;

coh_sim = sim_dat.dataMatSim(:, 1);     % Coherence of simulated trials
choice_sim = sim_dat.dataMatSim(:, 2);  % Choice of simulated trials
rt_sim = sim_dat.dataMatSim(:, 3);      % Reaction time of simulated trials

t = round(medA.leverage(1).t, 3);       % Time points used for original mediation analysis
dt = round(t(2) - t(1), 3);             % dt of original mediation analysis
% plotT = find(t >=0.2 & t <=0.5);
plot_tpts = t(find(t >=start_t & t <= end_t));

t_sim = round(sim_dat.tSim + 0.2, 3); % Time point 0 in tSim equalt 0.2s after motion onset
dt_sim = round(t_sim(2) - t_sim(1), 3);

useTPT = [];
for i = 1 : length(plot_tpts)
    useTPT(i) = find(t_sim == plot_tpts(i));
end
t_simPlot = t_sim(useTPT);

trls = 1 : 5; colrs = {'r', 'b', 'g', 'k', 'c', 'y'};
figure; hold on
for tr = trls
    plot(t_sim, sim_dat.diffusionPaths(tr, :)', 'Color', colrs{tr})
    line(rt_sim(tr)*[1 1], [-0.2 1.2], 'Color', colrs{tr})
end
title('t_ sim = tSim + 0.2')
xlabel('t_ sim')
ylabel('diffusionPaths')


%% Compute mediation on a subsample of trials with added noise to account for limited number of neurons
% 1000 permutations of selecting a subset of trials from the entire simulated data

I_orig = find(abs(coh_sim) < max_coh);      % include coherences 0%, 3.2% and 6.4%
t_mediation = find(t_sim(useTPT) == 0.55);
convSim = conv2(1, filtW50b, sim_dat.diffusionPaths, 'same'); % convolve with 50ms boxcar filter
sProj_sim_orig = convSim(:, useTPT);  
sProj_mediator_sim_orig = sProj_sim_orig(:,t_mediation);

numPerm = 1000;

% Scaling of noise for each neuron and time point
std_noise = 1.03;       % Bound height B = 1.03 for monkey M

% Settings for corr_with_RT_choice
do_plot = 0; % 0: Don't plot figures within function
flags.norm_to_se = 1;
noise_type = 'gauss';
noise_type = 'qrand'; rmin = 0.02; rmax = 0.08; nocheck = 0; std_noise = 1; 
noise_type = 'qrand'; rmin = 0.08; rmax = 0.09; nocheck = 0; std_noise = 1; 
noise_type = 'qrand'; rmin = 0; rmax = 0.09; nocheck = 0; std_noise = 1; 
% noise_type = 'qrand'; rmin = -2*10e-5; rmax = 7*10e-5; nocheck = 0; std_noise = 1; 
% too low
warning('off','all')

clear rho_RT_sess rho_RT_partial_sess rho_choice_sess rho_choice_partial_sess
clear rho_RT_sess_Sig rho_RT_partial_sess_Sig rho_choice_sess_Sig rho_choice_partial_sess_Sig
for sess = 7 : 8
    
    disp(num2str(sess))
    
    beh = par_sess{sess}.beh;   % structure of behavioral indicators of this session
    
    numTrials = length(find(abs(beh.sig_coh) < 0.1));   % number of trials with coh < 0.1 in this session
    numN = length(par_sess{sess}.unit_Tin);             % number of Tin_con neurons in this session
    
    for perm = 1 : numPerm
        % Select numTrials random trials to include in the analysis for this session
        if rem(perm, 50) == 0; disp(['Sess ' num2str(sess) ', perm ' num2str(perm)]); end
        this_perm = randperm(length(I_orig));
        this_perm = this_perm(1, 1 : numTrials);
        I = I_orig(this_perm);
        
        I_sess{sess, perm} = I;
        
        % Simulated 'true' DV (trials x time points)
        signal = sProj_sim_orig(I, :);
        
        % Generate numN noisy neurons in this session
        switch noise_type
            case 'gauss'
                % Add random noise to each neuron, trial and time window
                noise = std_noise * randn(numN, numTrials, size(signal, 2));
                sProj_sim = signal + squeeze(nanmean(noise, 1));
            case 'qrand'
                [C, Q] = qrancorrelmtx(numN,rmin,rmax,nocheck);
                
                % Q* signal
                sig_in = repmat(signal, 1, 1, numN);
                M_sig = permute(sig_in, [3 1 2]);
                QM_sig = squeeze(mmx('mult', Q, M_sig));
                sProj_simSig = squeeze(nanmean(QM_sig, 1));
                
                % Q *(signal + noise)
                noise = std_noise * randn(numN, numTrials, size(signal, 2));
                M_noise = M_sig + noise;
                QM_noise = squeeze(mmx('mult', Q, M_noise));
                sProj_sim = squeeze(nanmean(QM_noise, 1));
                
                % signal + Q * noise
                noise = std_noise * randn(numN, numTrials, size(signal, 2));
                M_sig_Qnoise = M_sig + squeeze(mmx('mult', Q, noise));
                sProj_simNoise = squeeze(nanmean(M_sig_Qnoise, 1));
                
                for trl = 1 : numTrials
                    C = cov(squeeze(QM_sig(:, trl, :))', 'partialrows');
                    cov_sess{sess}(:, trl, :, 1, perm) = C;
                    C = cov(squeeze(QM_noise(:, trl, :))', 'partialrows');
                    cov_sess{sess}(:, trl, :, 2, perm) = C;
                    C = cov(squeeze(M_sig_Qnoise(:, trl, :))', 'partialrows');
                    cov_sess{sess}(:, trl, :, 3, perm) = C;
                end
                
        end
        % Signal + noise
        sProj_mediator_sim = sProj_sim(:, t_mediation);
        
        [~,out_sim] = corr_with_RT_choice(choice_sim(I)==0, rt_sim(I), coh_sim(I), plot_tpts, ...
            sProj_sim, sProj_mediator_sim, minRT, maxRT, start_t, end_t, do_plot, flags);
        
        rho_RT_sess(sess, perm, :) = out_sim.rho_RT;
        rho_RT_partial_sess(sess, perm, :) = out_sim.rho_RT_partial;
        rho_choice_sess(sess, perm, :) = out_sim.rho_choice;
        rho_choice_partial_sess(sess, perm, :) = out_sim.rho_choice_partial;
        
        % Signal 
        sProj_mediator_simSig = sProj_simSig(:, t_mediation);
        
        [~,out_sim] = corr_with_RT_choice(choice_sim(I)==0, rt_sim(I), coh_sim(I), plot_tpts, ...
            sProj_simSig, sProj_mediator_simSig, minRT, maxRT, start_t, end_t, do_plot, flags);
        
        rho_RT_sess_Sig(sess, perm, :) = out_sim.rho_RT;
        rho_RT_partial_sess_Sig(sess, perm, :) = out_sim.rho_RT_partial;
        rho_choice_sess_Sig(sess, perm, :) = out_sim.rho_choice;
        rho_choice_partial_sess_Sig(sess, perm, :) = out_sim.rho_choice_partial;
        
        % Signal Q Noise
        sProj_mediator_simNoise = sProj_simNoise(:, t_mediation);
        
        [~,out_sim] = corr_with_RT_choice(choice_sim(I)==0, rt_sim(I), coh_sim(I), plot_tpts, ...
            sProj_simNoise, sProj_mediator_simNoise, minRT, maxRT, start_t, end_t, do_plot, flags);
        
        rho_RT_sess_Noise(sess, perm, :) = out_sim.rho_RT;
        rho_RT_partial_sess_Noise(sess, perm, :) = out_sim.rho_RT_partial;
        rho_choice_sess_Noise(sess, perm, :) = out_sim.rho_choice;
        rho_choice_partial_sess_Noise(sess, perm, :) = out_sim.rho_choice_partial;
        
    end
end
cov_sess_sim = cov_sess;
save(fullfile(saveLoc, 'mediation_sim_Q_0-009'), 'cov_sess_sim','rmin', 'rmax', 'nocheck',...
    'rho_RT_sess', 'rho_RT_partial_sess', 'rho_choice_sess', 'rho_choice_partial_sess',...
    'rho_RT_sess_Sig', 'rho_RT_partial_sess_Sig', 'rho_choice_sess_Sig', 'rho_choice_partial_sess_Sig',...
    'rho_RT_sess_Noise', 'rho_RT_partial_sess_Noise', 'rho_choice_sess_Noise', 'rho_choice_partial_sess_Noise', 'all_x', '-v7.3')


% Covariance between neurons in 3 signal + noise combinations
histbins = [0.13 : 0.001 : 0.15];
titles = {'Q*signal', 'Q* (signal + noise)', 'signal + Q*noise'};
figure; hold on
m = 3;
all_x = [];
for sess = 1 : 8
    disp(num2str(sess))
    A = squeeze(nanmean(nanmean(cov_sess{sess}(:, :, :, m, :), 5), 2));
    
    r = cov2r(A);
    LTril_cov = tril(r, -1);
    %         LTril_cov = tril(A, -1);
    LTril_cov(LTril_cov==0) = nan;
    
    x = LTril_cov(:);
    all_x = [all_x x'];
    
end
title(titles{m})
hist(x(~isnan(x)), histbins) %, [0.06 : 0.001 : 0.09])
xlim(histbins([1 end]))
xlabel('covariance');
if m == 1; ylabel('# simulated neuron pairs'); end
suptitle('distribution of pairwise correlations (simulation)')

histbins = []; 
titles = {'Q*signal', 'Q* (signal + noise)', 'signal + Q*noise'};
figure; hold on
for m = 1 : 3
    all_x = [];
    for sess = 1 : 8
        disp(num2str(sess))
        A = squeeze(nanmean(nanmean(cov_sess{sess}(:, :, :, m, :), 5), 2));
        
        r = cov2r(A);
        LTril_cov = tril(r, -1);
%         LTril_cov = tril(A, -1);
        LTril_cov(LTril_cov==0) = nan;
        
        x = LTril_cov(:);
        all_x = [all_x x'];
        
    end
    subplot(1, 3, m); hold on
    title(titles{m})
    hist(x(~isnan(x))) %, [0.06 : 0.001 : 0.09])
%     xlim([0.06 0.09])
    xlabel('covariance'); 
    if m == 1; ylabel('# simulated neuron pairs'); end
    
end
suptitle('distribution of pairwise correlations (simulation)')




figure; hold on
subplot(1, 2, 1); hold on
imagesc(Q1, [0 0.05])
subplot(1, 2, 2); hold on
imagesc(Q2, [0 0.05])

figure; hold on
trls = 1:5; colms = length(trls); ylimi = [-2.5 2.5];
for trl = 1 : colms
    ttt = trls(trl);
   subplot(4, colms, trl); hold on
   plot(squeeze(nanmean(QM_sig1(:, ttt, :), 1)), 'k', 'Linewidth', 2)
   plot(squeeze(QM_sig1(:, ttt, :))')
   ylim(ylimi)
   subplot(4, colms, trl+colms); hold on
   plot(squeeze(nanmean(QM_sig2(:, ttt, :), 1)), 'k', 'Linewidth', 2)
   plot(squeeze(QM_sig2(:, ttt, :))')
   
   ylim(ylimi)
   subplot(4, colms, trl+2*colms); hold on
   plot(squeeze(QM_noise1(:, ttt, :))')
   plot(squeeze(nanmean(QM_noise1(:, ttt, :), 1)), 'k', 'Linewidth', 2)
   ylim(ylimi)
   subplot(4, colms, trl+3*colms); hold on
   plot(squeeze(QM_noise2(:, ttt, :))')
   plot(squeeze(nanmean(QM_noise2(:, ttt, :), 1)), 'k', 'Linewidth', 2)
   ylim(ylimi)
end


% C_sig = cov(squeeze(QM_sig(:, 1, :))');
% C_noise = cov(squeeze(QM_noise(:, 1, :))');
% 
% figure; hold on
% subplot(1, 2, 1); imagesc(C_sig, [-0.2 0.2])
% subplot(1, 2, 2); imagesc(C_noise, [-0.2 0.2])
% 
% 
% figure; hold on
% for trl = 1 : 4
%    subplot(2, 2, trl); hold on 
% %    plot(squeeze(mmx('mult', Q, M_noise(:, trl, :)))', 'r')
% %    plot(squeeze(mmx('mult', Q, M_sig(:, trl, :)))', 'b')
%    
%    plot(signal(trl, :), 'k', 'LineWidth', 2)
%    plot(sProj_sim(trl, :), 'm', 'LineWidth', 2)
%    plot(sProj_sim2(trl, :), 'b', 'LineWidth', 2)
%    plot(sProj_simSig(trl, :), 'c', 'LineWidth', 2)
%    
% end


plot_med = find(t_simPlot<=0.50);
figure; hold on
subplot(1, 2, 1); hold on
plot(t_simPlot, squeeze(rho_choice_sess(sess, perm, :)), 'r')
plot(t_simPlot, squeeze(rho_choice_sess_Sig(sess, perm, :)), 'b')
plot(t_simPlot(plot_med), squeeze(rho_choice_partial_sess(sess, perm, plot_med)), 'r--')
plot(t_simPlot(plot_med), squeeze(rho_choice_partial_sess_Sig(sess, perm, plot_med)), 'b--')
plot(t_simPlot, squeeze(rho_choice_sess_Noise(sess, perm, :)), 'g')
plot(t_simPlot(plot_med), squeeze(rho_choice_partial_sess_Noise(sess, perm, plot_med)), 'g--')
subplot(1, 2, 2); hold on
plot(t_simPlot, squeeze(rho_RT_sess(sess, perm, :)), 'r')
plot(t_simPlot, squeeze(rho_RT_sess_Sig(sess, perm, :)), 'b')
plot(t_simPlot(plot_med), squeeze(rho_RT_partial_sess(sess, perm, plot_med)), 'r--')
plot(t_simPlot(plot_med), squeeze(rho_RT_partial_sess_Sig(sess, perm, plot_med)), 'b--')
plot(t_simPlot, squeeze(rho_RT_sess_Noise(sess, perm, :)), 'g')
plot(t_simPlot(plot_med), squeeze(rho_RT_partial_sess_Noise(sess, perm, plot_med)), 'g--')

% save(fullfile(saveLoc, 'mediation_sim_noise_perm'), 'I_sess', 'rho_RT_sess', 'rho_RT_partial_sess', 'rho_choice_sess', 'rho_choice_partial_sess',...
%     'minRT', 'maxRT', 'start_t', 'end_t')


% Plot
lwR = 2;
lwM = 1;
col_self = [48 136 193]/255;
col_Tin = [237 190 4]/255;
xLims = [0.2 0.56]; 
yLimCho = [-2 10];
yLimRT = [-0.4 0.1];
sFont = 6;
plotT = plot_med; 

figure('Position', [200 10 700 250]); hold on
% Choice =============================================================
% Leverage on choice: neuron classes
subplot(1, 2, 1); hold on
% Raw choice
m = squeeze(nanmean(nanmean(rho_choice_sess, 1), 2));
se = squeeze(nanstd(nanmean(rho_choice_sess, 1), 1));
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(t_mediation), m(t_mediation), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(t_mediation)*[1 1], m(t_mediation) + se(t_mediation) * [-1 1], 'Color', plotCol)
% Mediated by self
m = squeeze(nanmean(nanmean(rho_choice_partial_sess, 1), 2));
se = squeeze(nanstd(nanmean(rho_choice_partial_sess, 1), 1));
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
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
m = squeeze(nanmean(nanmean(rho_RT_sess, 1), 2));
se = squeeze(nanstd(nanmean(rho_RT_sess, 1), 1));
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(t_mediation), m(t_mediation), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(t_mediation)*[1 1], m(t_mediation) + se(t_mediation) * [-1 1], 'Color', plotCol)
% Mediated by self
m = squeeze(nanmean(nanmean(rho_RT_partial_sess, 1), 2));
se = squeeze(nanstd(nanmean(rho_RT_partial_sess, 1), 1));
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwM;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
% Formatting
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim([-0.6 0.1])
xlabel('Time from motion onset [s]');
ylabel('Correlation with RT (r)')
set(gca,'TickDir','out', 'FontSize', sFont);
suptitle('Q * (sinal + noise)')




%%

figure('Position', [200 10 700 250]); hold on
% Choice =============================================================
% Leverage on choice: neuron classes
subplot(1, 2, 1); hold on
% Raw choice
m = squeeze(nanmean(nanmean(rho_choice_sess_Sig, 1), 2));
se = squeeze(nanstd(nanmean(rho_choice_sess_Sig, 1), 1));
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(t_mediation), m(t_mediation), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(t_mediation)*[1 1], m(t_mediation) + se(t_mediation) * [-1 1], 'Color', plotCol)
% Mediated by self
m = squeeze(nanmean(nanmean(rho_choice_partial_sess_Sig, 1), 2));
se = squeeze(nanstd(nanmean(rho_choice_partial_sess_Sig, 1), 1));
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
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
m = squeeze(nanmean(nanmean(rho_RT_sess_Sig, 1), 2));
se = squeeze(nanstd(nanmean(rho_RT_sess_Sig, 1), 1));
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(t_mediation), m(t_mediation), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(t_mediation)*[1 1], m(t_mediation) + se(t_mediation) * [-1 1], 'Color', plotCol)
% Mediated by self
m = squeeze(nanmean(nanmean(rho_RT_partial_sess_Sig, 1), 2));
se = squeeze(nanstd(nanmean(rho_RT_partial_sess_Sig, 1), 1));
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwM;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
% Formatting
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim([-0.6 0.1])
xlabel('Time from motion onset [s]');
ylabel('Correlation with RT (r)')
set(gca,'TickDir','out', 'FontSize', sFont);
suptitle('Q * sinal')


%% Signal + Q * noise

figure('Position', [200 10 700 250]); hold on
% Choice =============================================================
% Leverage on choice: neuron classes
subplot(1, 2, 1); hold on
% Raw choice
m = squeeze(nanmean(nanmean(rho_choice_sess_Noise, 1), 2));
se = squeeze(nanstd(nanmean(rho_choice_sess_Noise, 1), 1));
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(t_mediation), m(t_mediation), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(t_mediation)*[1 1], m(t_mediation) + se(t_mediation) * [-1 1], 'Color', plotCol)
% Mediated by self
m = squeeze(nanmean(nanmean(rho_choice_partial_sess_Noise, 1), 2));
se = squeeze(nanstd(nanmean(rho_choice_partial_sess_Noise, 1), 1));
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
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
m = squeeze(nanmean(nanmean(rho_RT_sess_Noise, 1), 2));
se = squeeze(nanstd(nanmean(rho_RT_sess_Noise, 1), 1));
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(t_mediation), m(t_mediation), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(t_mediation)*[1 1], m(t_mediation) + se(t_mediation) * [-1 1], 'Color', plotCol)
% Mediated by self
m = squeeze(nanmean(nanmean(rho_RT_partial_sess_Noise, 1), 2));
se = squeeze(nanstd(nanmean(rho_RT_partial_sess_Noise, 1), 1));
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwM;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
% Formatting
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim([-0.6 0.1])
xlabel('Time from motion onset [s]');
ylabel('Correlation with RT (r)')
set(gca,'TickDir','out', 'FontSize', sFont);
suptitle('signal + Q * noise')

%%




figure; hold on
subplot(1, 2, 1); hold on
plot(squeeze(nanmean(rho_RT_partial_sess_Sig(1:5, 1:5, :), 1))')
subplot(1, 2, 2); hold on
plot(squeeze(nanmean(rho_RT_partial_sess(1:5, 1:5, :), 1))')







