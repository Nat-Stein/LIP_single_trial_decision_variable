addpath('../generic/');
addpath(genpath('../az_matlab_files/'));

%%

ndatasets = 8;

names = {'all'};
vlambda = [0.0001,0.001,0.005,0.01,0.02,0.05,0.1];
nlambda = length(vlambda);
AUC = nan(ndatasets, length(names),nlambda);
PROP_zeros = AUC;
clear vB

for idataset = 1:ndatasets

    dt_ms = 25; % 100
    %     dt_ms = 50;
    dataset = idataset;
    dat = get_data(dt_ms,[],dataset,'dt_rel_RT',0.05);

    struct2vars(dat);

    idx_include = choice==0;
    %     idx_include = choice==1;

    dt = t(2)-t(1);

    R = [];
    D = [];
    trial_id = [];
    for i=1:length(t)-1
        r = zeros(size(RT));
        r(RT>=t(i) & RT<(t(i) + 0.15)) = 1;
        R = [R; r(idx_include)];

        D = cat(1,D,squeeze(H(:,i,idx_include))');

        trial_id = [trial_id; find(idx_include)];
    end

    % remove nans
    I = all(isnan(D),2);
    R = R(~I);
    D = D(~I,:);
    trial_id = trial_id(~I);

    %% regre

    allneurons = [1:dat.nneurons]';
    idx_neurons = {allneurons};

    J = iseven(trial_id);
    %     W = zeros(dat.nneurons,length(idx_neurons));
    for k=1:length(idx_neurons)

        if ~isempty(idx_neurons{k})

            for ilambda = 1:nlambda

                [B,FitInfo] = lassoglm(D(J,idx_neurons{k}),R(J),'binomial','link','logit','lambda',vlambda(ilambda));
                B0 = FitInfo.Intercept;

                coef = [B0; B];

                % pred
                yhat = glmval(coef,D(~J,idx_neurons{k}),'logit');

                %% calc ROC
                scores = yhat;
                labels = R(~J);
                posclass = 1;
                [~,~,~,auc] = perfcurve(labels,scores,posclass);
                AUC(idataset,k,ilambda) = auc;
                PROP_zeros(idataset,k,ilambda) = mean(B==0);
                vB{idataset,ilambda,k} = B;
            end
            %%
            % %             W(idx_neurons{k},k) = B;
            %             Weights{k} = B; % regression weights
        end

    end
    %%

    %     savefilename = fullfile('./data_by_neuronal_group',dataset);
    %     save(savefilename,'W','names');


end

%%
p = publish_plot(3,1);

p.next();
plot(vlambda,nanmean(squeeze(AUC)),'.-');
set(gca,'xscale','log','xtick',vlambda,'xticklabel',vlambda,'xminortick','off');
xlabel('Lambda');
ylabel('AUC on left-out trials');

p.next();
plot(vlambda,nanmean(squeeze(PROP_zeros)),'.-');
set(gca,'xscale','log','xtick',vlambda,'xticklabel',vlambda,'xminortick','off');
xlabel('Lambda');
ylabel('Proportion of neurons with zero weight');

p.next();
ind = find(vlambda==0.01);
for i=1:nlambda
    for j=1:ndatasets
        coSim(j,i) = 1 - pdist2(vB{j,i}',vB{j,ind}','cosine');
    end
end
plot(vlambda,nanmean(coSim),'.-');
set(gca,'xscale','log','xtick',vlambda,'xticklabel',vlambda,'xminortick','off');
xlabel('Lambda');
ylabel('Cosine similarity to lambda = 0.01');
ylim([0,1]);
% a = vB{1,4}';
% b = vB{1,2}';
% cosSim = dot(a,b)/(norm(a)*norm(b)); 

p.format('FontSize',12);
p.append_to_pdf('fig_crossvalid',1,1);



