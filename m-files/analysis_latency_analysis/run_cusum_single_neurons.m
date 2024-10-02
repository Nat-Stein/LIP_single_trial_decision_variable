% addpath('../generic/');
% addpath('../analysis_reclassif/');

%%

n_datasets = 8;
fig_made = 0;
for idataset = 1 : n_datasets

%     addpath('../analysis_explore/');
    dt_ms = 25;
    par = par_sess{idataset};

    D = get_data(dt_ms,[],idataset, par);

    %     struct2vars(D);
    filt = ones(size(D.RT))==1;

    %%
    latency_cusum = nan(D.nneurons,1);
    AUC = nan(D.nneurons,1);

    for id = 1 : D.nneurons


        tr_ind = abs(D.coh)>0.1 & D.correct==1 & D.RT>0.45;
        h = squeeze(D.H(id,:,tr_ind));

        % calc ROC
        % direction selctive
        tind = D.t>0.1 & D.t<0.4;
        scores = nanmean(h(tind,:),1)';
        labels = D.coh(tr_ind)>0;
        posclass = 0;
        [X,Y,T,auc,OPTROCPT,SUBY] = perfcurve(labels,scores,posclass);

        AUC(id) = max(auc, 1-auc);

        h = squeeze(nanmean(D.H(id,:,D.coh<0 & D.correct==1),3)) - ...
            squeeze(nanmean(D.H(id,:,D.coh>0 & D.correct==1),3));

        % cusum
        tind  = find(D.t>=0 & D.t<=0.5); % for cusum

        tt = D.t(tind);
        ch = cumsum(h(tind));

        ind = length(tind);

        % is there a significant modulation by motion coh? If so, calc latency

        % calc latency
        if ~all(isnan(h(tind)))
            

            [lat,i0,a,b,d,m] = cusum.cusum(tt(1:ind),ch(1:ind),do_plot,0);

            latency_cusum(id) = lat;
           
        end
    end

    %% save
    minCells = D.minCells;
    nneurons = D.nneurons;
    unitIdxLIP_Tin = D.unitIdxLIP_Tin;
    unitIdxLIP_Tout = D.unitIdxLIP_Tout;
    here= pwd;
    saveLoc_latency = fullfile(here, 'analysis_latency_analysis');
    save(fullfile(saveLoc_latency, 'data', D.dataset),'AUC','latency_cusum','minCells','unitIdxLIP_Tin','unitIdxLIP_Tout','nneurons');
    
end