% Set up format of simulation data to match mediation analysis

medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));
% sim_dat = load(fullfile(saveLoc, 'race_model_simulation'));


% use race model with lower bound to simulate traces
% [S,P]  = simRace2DDM(t, varargin)
numTr = 60000;

t_sim = (0:0.001:5);
cohs_sim = [cohs 0]; % simulate two sets of 0% coherence

clear sim_dat
sim_dat.tSim = t_sim;
ct = 0;
for coh = cohs_sim
    disp(num2str(coh))
    for n = 1 : numTr / length(cohs_sim)
        ct = ct + 1;
        [X,P] = simRace2DDM(t_sim, 'coh', coh);
        simRace2DDM(t_sim,'seed',P.seed);
        
        sim_dat.dataMatSim(ct, :) = [coh X.choice X.decTime];
        sim_dat.diffusionPaths(ct, :) = X.dvTin;
        
    end
end
sim_dat.dataMatSim(:, 2) = 1 - sim_dat.dataMatSim(:, 2);
% save(fullfile(saveLoc, 'race_model_simulation_mns'), 'sim_dat', 'P', 't', '-v7.3');

%% normalize simulated signal
for trl = 1 : numTr
    [~, index] = min(abs(t_sim-sim_dat.dataMatSim(trl, 3)));
    atRT(trl) = sim_dat.diffusionPaths(trl, index);
end
mm = [];
for cho = 0 : 1
    trls = sim_dat.dataMatSim(:, 2) == cho;
    sl = nanmean(sim_dat.diffusionPaths(trls, :), 1);
    rl = nanmean(atRT(trls));
    mm = [ns sl rl];
end

maxM = max(mm);
minM = min(mm);
exc = maxM - minM;

sim_dat.diffusionPaths_norm = (sim_dat.diffusionPaths - minM) ./ exc;

save(fullfile(saveLoc, 'race_model_simulation_mns'), 'sim_dat', 'P', '-v7.3');


%% Plotting without noise

start_t = sim_dat.tSim(1);
end_t = 0.6;    % was 0.55 before

max_coh = 0.1;
minRT = 0.67;
maxRT = 2;

coh_sim = sim_dat.dataMatSim(:, 1);
choice_sim = sim_dat.dataMatSim(:, 2);
rt_sim = sim_dat.dataMatSim(:, 3);

t = round(medA.leverage(1).t, 3);
dt = round(t(2) - t(1), 3); 
plotT = find(t >=0.2 & t <=0.5);
plot_tpts = t(find(t >=start_t & t <= end_t));

t_sim = round(sim_dat.tSim, 3);
dt_sim = round(t_sim(2) - t_sim(1), 3);

for i = 1 : length(plot_tpts)
    useTPT(i) = find(t_sim == plot_tpts(i));
end

%%  -----

I = find(abs(coh_sim) < 0.1);

sProj_sim = sim_dat.diffusionPaths_norm(:, useTPT);  % convolve with correct filter
sProj_mediator_sim = sProj_sim(:,end);

% Compute mediation
[~,out_sim] = corr_with_RT_choice(choice_sim(I)==0, rt_sim(I), coh_sim(I), plot_tpts, ...
    sProj_sim(I,:), sProj_mediator_sim(I), minRT, maxRT, start_t, end_t, do_plot,flags);

% Plot mediation of simulated data
plotT = find(out_sim.tt >=0.2 & out_sim.tt <=0.5);

figure('Position', [200 10 700 250]); hold on
% Choice =============================================================
% Leverage on choice: neuron classes
subplot(1, 2, 1); hold on
plotCol = [0 0 0];
plot(out_sim.tt(plotT), out_sim.rho_choice(plotT), 'Color', plotCol, 'LineWidth', 2)
plot(out_sim.tt(end), out_sim.rho_choice(end), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
plotCol = col_self;
plot(out_sim.tt(plotT), out_sim.rho_choice_partial(plotT), 'Color', plotCol, 'LineWidth', 2)
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim([-10 30])
xlabel('Time from motion onset [s]');
ylabel('Leverage on choice (slope beta)')
set(gca,'TickDir','out', 'FontSize', sFont);
% RT ==================================================================
% Correlation with RT
subplot(1, 2, 2); hold on
plotCol = [0 0 0];
plot(out_sim.tt(plotT), out_sim.rho_RT(plotT), 'Color', plotCol, 'LineWidth', 2)
plot(out_sim.tt(end), out_sim.rho_RT(end), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
plotCol = col_self;
plot(out_sim.tt(plotT), out_sim.rho_RT_partial(plotT), 'Color', plotCol, 'LineWidth', 2)
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim([-0.4 0.1])
xlabel('Time from motion onset [s]');
ylabel('Correlation with RT (r)')
set(gca,'TickDir','out', 'FontSize', sFont);



%% Compute mediation on a subsample of trials with added noise to account for limited number of neurons

I_orig = find(abs(coh_sim) < 0.1);
t_mediation = find(sim_dat.tSim(useTPT) == 0.55);
sProj_sim_orig = sim_dat.diffusionPaths_norm(:, useTPT);  % convolve with correct filter
sProj_mediator_sim_orig = sProj_sim_orig(:,t_mediation);

% Next different noise method: Q = sqrtm([1 P.r; P.r 1]);

do_plot = 0; 
flags.norm_to_se = 1;
% create 8 different sample sessions with trial numbers that match the data
std_noise = 0.9; %P.B; %1;
noise_source = 'weakC'; % 'poiss', 'gauss', 'weakC'
% Settings for 'weakC'
rmin = 0.02; rmax = 0.08; % too low
nocheck = 0;

clear sim_leverage
for sess = 1 : 8
    
    disp(num2str(sess))
    
    beh = par_sess{sess}.beh;   % structure of behavioral indicators of this session
    
    numTrials = length(find(abs(beh.sig_coh) < 0.1));   % number of trials with coh < 0.1 in this session
    numN = length(par_sess{sess}.unit_Tin);
    
    [C, Q] = qrancorrelmtx(numN,rmin,rmax,nocheck);
    
    this_perm = randperm(length(I_orig));
    this_perm = this_perm(1, 1 : numTrials);
    I = I_orig(this_perm);      % the simulated trials being used as a comparison for this session
    
    I_sess{sess} = I;
    
    signal = sProj_sim_orig(I, :);
    
    switch noise_source
        case 'poiss'
            % Poisson noise neurons
            randMat = std_noise * randn(numN, numTrials, size(signal, 2));
            for trl = 1 : numTrials
                clear spks
                for n = 1 : numN
                    spks(n, :) = signal(trl, :) > squeeze(randMat(n, trl, :))';
                end
                ss = mean(spks, 1); ss(isnan(signal(trl, :))) = nan;
                signal_meanPoiss(trl, :) = ss;
            end
            % noise = poissrnd(1, numN, numTrials, size(signal, 2));
            sProj_sim = signal_meanPoiss;
        case 'gauss'
            noise = std_noise * randn(numN, numTrials, size(signal, 2));
            sProj_sim = signal + squeeze(nanmean(noise, 1));
        case 'weakC'
            % Figure out which matrix has to be a transposed
            noise = std_noise * randn(numN, numTrials, size(signal, 2));
            
            sig_in = repmat(signal, 1, 1, numN);
            
                    M = permute(sig_in, [3 1 2]) + noise;
%             M_2 = noise;
            for trl = 1 : numTrials
                            M1 = squeeze(M(:, trl, :))';
                            M2 = M1*Q;
                            sProj_sim(trl, :) = nanmean(M2, 2);
                
                
%                 M1_2 = squeeze(M_2(:, trl, :))';
%                 M2_2 = M1_2*Q;
%                 sProj_sim(trl, :) = nanmean(M2_2, 2) + signal(trl, :)';
                
            end
            %         figure; hold on;
            %         plot(squeeze(signal(trl, :)))
            %         plot(squeeze(nanmean(M(:, trl, :), 1))')
            %         plot(squeeze(nanmean(M2(:, :), 2))')
            %         plot(nanmean(M2_2, 2) + signal(trl, :)')
        otherwise
            disp('noise case does not exist')
    end
    
    sProj_mediator_sim = sProj_sim(:, t_mediation);
    
    [~,out_sim] = corr_with_RT_choice(choice_sim(I)==0, rt_sim(I), coh_sim(I), plot_tpts, ...
        sProj_sim, sProj_mediator_sim, minRT, maxRT, start_t, end_t, do_plot,flags);
    
    sim_leverage.rho_RT_sess(sess, :) = out_sim.rho_RT;
    sim_leverage.rho_RT_partial_sess(sess, :) = out_sim.rho_RT_partial;
    sim_leverage.rho_choice_sess(sess, :) = out_sim.rho_choice;
    sim_leverage.rho_choice_partial_sess(sess, :) = out_sim.rho_choice_partial;
end

% save(fullfile(saveLoc, 'mediation_sim_noise'), 'I_sess', 'rho_RT_sess', 'rho_RT_partial_sess', 'rho_choice_sess', 'rho_choice_partial_sess',...
%     'minRT', 'maxRT', 'start_t', 'end_t')

% save(fullfile(saveLoc, 'mediation_sim'), 'I_sess', 'rho_RT_sess', 'rho_RT_partial_sess', 'rho_choice_sess', 'rho_choice_partial_sess',...
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

figure('Position', [200 10 700 250]); hold on
% Choice =============================================================
% Leverage on choice: neuron classes
subplot(1, 2, 1); hold on
% Raw choice
m = nanmean(sim_leverage.rho_choice_sess, 1);
se = nanstd(sim_leverage.rho_choice_sess, 1)/sqrt(8);
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(end), m(end), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(end)*[1 1], m(end) + se(end) * [-1 1], 'Color', plotCol)
% Mediated by self
m = nanmean(sim_leverage.rho_choice_partial_sess, 1);
se = nanstd(sim_leverage.rho_choice_partial_sess, 1)/sqrt(8);
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwM;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
% Formatting
line(xLims, [0 0], 'Color', 'k')
xlim(xLims)
ylim(yLimCho)
xlabel('Time from motion onset [s]');
ylabel('Leverage on choice (slope beta)')
set(gca,'TickDir','out', 'FontSize', sFont);
% RT ==================================================================
% Correlation with RT
subplot(1, 2, 2); hold on
plotCol = [0 0 0];
% Raw choice
m = nanmean(sim_leverage.rho_RT_sess, 1);
se = nanstd(sim_leverage.rho_RT_sess, 1)/sqrt(8);
plotCol = 'k';
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
set(s.edge,'LineWidth',0.5,'LineStyle','-', 'Color', plotCol)
s.mainLine.LineWidth = lwR;
s.patch.FaceColor = plotCol;
s.mainLine.Color = plotCol;
plot(out_sim.tt(end), m(end), 'o', 'Color', plotCol, 'MarkerFaceColor', plotCol)
line(out_sim.tt(end)*[1 1], m(end) + se(end) * [-1 1], 'Color', plotCol)
% Mediated by self
m = nanmean(sim_leverage.rho_RT_partial_sess, 1);
se = nanstd(sim_leverage.rho_RT_partial_sess, 1)/sqrt(8);
plotCol = col_self;
s = shadedErrorBar(out_sim.tt(plotT), m(plotT), se(plotT));
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






