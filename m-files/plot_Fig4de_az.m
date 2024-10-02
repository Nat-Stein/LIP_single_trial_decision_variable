
% Load data

load(fullfile(saveLoc, 'cosine_sim'), 'cosine_sim')
% load(fullfile(saveLoc, 'resCorr_wiTrial'), 'rz_sess');

load(fullfile(saveLoc, 'single_trial_corr'))
rho = cat(3,C.rho);
rz_sess = nanmean(rho,3);

% load(fullfile(saveLoc, 'cosine_sim_randperm'), 'cosine_sim_randperm') % cosine_sim_randperm(sess, d1, d2, :)

% Plot d
load('E:\Matalb analyses\darkRedBlue')
figure('Position', [600 400 800 400]); hold on;
subplot(1, 2, 1); hold on; % title('Cosine similarity'); 
A = squeeze(nanmean(cosine_sim, 1)); A = tril(A); A = flip(A);
imagesc(A, [-1 1]);
set(gca,'xtick',1:length(dimLabelsCS),'xticklabel', dimLabelsCS,'tickdir','out');
set(gca,'ytick',1:length(dimLabelsCS),'yticklabel', flip(dimLabelsCS),'tickdir','out');
xtickangle(45)
colormap(darkRedBlue);
cb = colorbar;
cb.Location = 'northoutside';
cb.Ticks = [-1 1];
axis tight
axis square
% Plot e
subplot(1, 2, 2); hold on; % title('r')
A = rz_sess; % A = r_sess(useDims, useDims); 
A = tril(A); A = flip(A);
imagesc(A, [-1 1]);
set(gca,'xtick',1:length(dimLabelsCS),'xticklabel', dimLabelsCS,'tickdir','out');
set(gca,'ytick',1:length(dimLabelsCS),'yticklabel', {'', '', '', '', ''},'tickdir','out');
xtickangle(45)
colormap(darkRedBlue);
cb = colorbar;
cb.Location = 'northoutside';
cb.Ticks = [-1 1];
axis tight
axis square



