

load(fullfile(saveLoc, 'xMediationMatrices'))


% load('E:\Matalb analyses\xMediationMatrices_230522')

% indY = [6 5 4 3 7];
% indX = [7 3 4 5 6];
indX = [6 4 5 7 8];
indY = [8 7 5 4 6];

indX = [7 5 6 8 9];
indY = flip(indX);

figure('Position', [100  100  900  400]); hold on
% ind = 3:length(uni_s);
subplot(1, 2, 1); hold on; title('Mediation of leverage on choice')
imagesc(X(indY, indX), [0 1]);
set(gca,'xtick',1:length(uni_s),'xticklabel',labels(indX),'tickdir','out');
set(gca,'ytick',1:length(uni_s),'yticklabel',labels(indY));
colores = cbrewer('seq','YlOrRd',100);
colormap(colores);
colorbar
axis tight
axis square
xlabel('Mediator S^y');
ylabel('Signal S^x');

subplot(1, 2, 2); hold on; title('Mediation of correlation with RT')
imagesc(Xrt(indY, indX), [0 1]);
set(gca,'xtick',1:length(uni_s),'xticklabel',labels(indX),'tickdir','out');
set(gca,'ytick',1:length(uni_s),'yticklabel',labels(indY));
colores = cbrewer('seq','YlOrRd',100);
colormap(colores);
colorbar
axis tight
axis square
xlabel('Mediator S^y');
ylabel('Signal S^x');


