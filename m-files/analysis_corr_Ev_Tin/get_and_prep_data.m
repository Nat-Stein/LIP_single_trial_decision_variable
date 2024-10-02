function [t,TIN, TOUT,EC,EI,CHOICE,RT,COH,SESSION] = get_and_prep_data()

% Adaptend from AZ 2022

TIN = [];
TOUT = [];
EC = [];
EI = [];
RT = [];
COH = [];
CHOICE = [];
SESSION = [];

% do_save_flag = 1;
n_datasets = 8;
global par_sess

for idataset = 1 : n_datasets

%     dt_ms = 10;
    dt_ms = 25;
    par = par_sess{idataset};

    D = get_data(dt_ms, [], idataset, par);
%     struct2vars(D);
    filt_trials = ones(size(D.choice))==1;

    %%
    %%
    I = D.minCells.DinRFc;
    K = D.minCells.DinRFi;
    J = D.unitIdxLIP_Tin;
    L = D.unitIdxLIP_Tout;

    idx = ismember(I,[J,L]);
    I(idx) = [];

    idx = ismember(K,[J,L]);
    K(idx) = [];
    
    sTin = squeeze(nanmean(D.H(J,:,:),1))';
    sEc = squeeze(nanmean(D.H(I,:,:),1))';
    sEi = squeeze(nanmean(D.H(K,:,:),1))';
    sTout = squeeze(nanmean(D.H(L,:,:),1))';

    %% calc residuals relative the signed coh
    ucoh = unique(D.coh);
    for i=1:length(ucoh)
        I = D.coh==ucoh(i) & filt_trials;
        sTin(I,:) = sTin(I,:) - nanmean(sTin(I,:));
        sEc(I,:) = sEc(I,:) - nanmean(sEc(I,:));
        sEi(I,:) = sEi(I,:) - nanmean(sEi(I,:));
        sTout(I,:) = sTout(I,:) - nanmean(sTout(I,:));
    end
    %% standarize
    sTin = nanzscore(sTin);
    sEc = nanzscore(sEc);
    sEi = nanzscore(sEi);
    sTout = nanzscore(sTout);

    TIN = [TIN; sTin];
    EC = [EC; sEc];
    EI = [EI; sEi];
    TOUT = [TOUT; sTout];

    RT =[RT; D.RT];
    COH = [COH; D.coh];
    CHOICE = [CHOICE; D.choice];
    SESSION = [SESSION; ones(size(D.RT))*idataset];

end

%% remove baseline
t = D.t;
tind = t>=-0.1 & t<=0;
EC = EC - nanmean(EC(:,tind),2);
EI = EI - nanmean(EI(:,tind),2);
TIN = TIN - nanmean(TIN(:,tind),2);
TOUT = TOUT - nanmean(TOUT(:,tind),2);

end
