addpath('../generic/');
addpath('../analysis_reclassif/');

%%

n_datasets = 8;
auc_thres = 0.6;

fields = {'DinRFc','DinRFi','unitIdxLIP_Tin','unitIdxLIP_Tout'};
inside = {'minCells','minCells','',''};
for i=1:length(fields)
    latency.(fields{i}) = [];
end

files = dir(fullfile(saveLoc_latency, 'data', '/*.mat'));
for idataset=1:length(files)

%     addpath('../analysis_explore/');
%     dt_ms = 100; % for speed, not used
%     D = get_data(dt_ms,[],idataset);
    
    D = load(fullfile(saveLoc_latency,files(idataset).name));
    
    latency.dataset(idataset).name = files(idataset).name;
    
    % group the latency of different neuron classes
    for i=1:length(fields)
        if isempty(inside{i})
            I = D.AUC>=auc_thres & ismember([1:D.nneurons]',D.(fields{i}));
        else
            I = D.AUC>=auc_thres & ismember([1:D.nneurons]',D.(inside{i}).(fields{i}));
        end
        latency.(fields{i}) = [latency.(fields{i}); D.latency_cusum(I)];
        
        latency.dataset(idataset).(fields{i}) = find(I);
    end

end

save(fullfile(saveLoc_latency, 'latency_estimates'), 'latency')

%% plot
t = linspace(0,0.5,1000);

figure; hold on
colores = movshon_colors(3);
lat_Tin = latency.unitIdxLIP_Tin;
lat_Din = [latency.DinRFc; latency.DinRFi];
plot(t,sum(t>lat_Tin)/length(lat_Tin),'color',colores(1,:));
hold all
plot(t,sum(t>lat_Din)/length(lat_Din),'color',colores(2,:));

m = mean(lat_Tin);
s = stderror(lat_Tin);
ha(1) = arrow([m,0.1],[m,0],'facecolor',colores(1,:),'edgecolor',colores(1,:));
he(1) = plot([m-s,m+s],[0,0],'color',colores(1,:),'LineWidth',2);

m = mean(lat_Din);
s = stderror(lat_Din);
ha(2) = arrow([m,0.1],[m,0],'facecolor',colores(2,:),'edgecolor',colores(2,:));
he(2) = plot([m-s,m+s],[0,0],'color',colores(2,:),'LineWidth',2);

legend('Target in','Dots in');

xlabel('Time from stim. onset [s]');
ylabel('Cumulative prob. of onset times across neurons')
set(ha,'LineWidth',2);
set(he,'LineWidth',3);

