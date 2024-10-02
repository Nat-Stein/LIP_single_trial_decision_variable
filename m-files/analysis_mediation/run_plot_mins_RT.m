load ../data/sigSL_allSessions_Ariel_230522.mat
S = sigSL_Ariel;
MinC_minus_I = S.MinC - S.MinI;

%%

RT = S.rt/1000;

ntr = length(S.choice);

% tind = S.t>=0.2 & S.t<=0.3;
% min_RT = 0.4;
% max_RT = 0.7;

tind = S.t>=0.1 & S.t<=0.4;
min_RT = 0.5;
max_RT = 2;

d = nansum(MinC_minus_I(:,tind),2);
% d = nansum(S.TinC(:,tind),2);
% d = nansum(S.MinC(:,tind),2);
% d = nansum(S.MinI(:,tind),2);
% d = nansum(S.whenC(:,tind),2);

% [~,~,u] = unique([S.choice, S.sig_coh],'rows');

% I = abs(S.sig_coh)<=0.128 & RT>min_RT;
% 
% conditions = [S.choice, S.sig_coh, S.session];
% [v(I)] = index_prctile_by_group(d(I),[0,50,100],conditions(I,:));

v = nan(ntr,1);
ucoh = nanunique(S.sig_coh);
for i=1:length(ucoh)
    for j=0:1
        for k=1:8 %sessions
            I = S.sig_coh==ucoh(i) & S.choice==j & RT>min_RT & RT<max_RT & S.session==k;
            if sum(I)>1
                mediana = nanmedian(d(I));
                v(I) = 1 + (d(I)>mediana);
            end
        end
    end
end

v2 = nan(ntr,1);
ucoh = nanunique(S.sig_coh);
for i=1:length(ucoh)
    for j=0:1
        for k=1:8 %sessions
            I = S.sig_coh==ucoh(i) & RT>min_RT & RT<max_RT & S.session==k;
            if sum(I)>1
                mediana = nanmedian(d(I));
                v2(I) = 1 + (d(I)>mediana);
            end
        end
    end
end

% debugging
% v2 = nan(ntr,1);
% ucoh = nanunique(S.sig_coh);
% for j=0:1
%     for k=1:8 %sessions
%         I = RT>min_RT & S.session==k;
%         if sum(I)>1
%             mediana = nanmedian(d(I));
%             v2(I) = 1 + (d(I)>mediana);
%         end
%     end
% end

%
p = publish_plot(2,1);
p.next();
J = true(size(S.choice));
curva_media(S.choice, 100*S.sig_coh,v2==1 & J,2);
hold all
curva_media(S.choice, 100*S.sig_coh,v2==2 & J,2);

ylabel('Prob. rightward choice');
h = legend('Low M_{contra} - M_{ipsi}','High M_{contra} - M_{ipsi}');

p.next();
curva_media(RT, 100*S.sig_coh, S.choice==0 & v==1 & J,2);
hold all
curva_media(RT, 100*S.sig_coh, S.choice==0 & v==2 & J,2);

xlabel('Motion strength (%)');
ylabel('Response time [s]');



same_xlim(p.h_ax);
p.format();
set(h,'location','best','fontsize',10,'box','off');

%%
p.append_to_pdf('fig_mins_RT',1,1);

%% stats
I = S.choice==0 & RT>min_RT;

depvar = RT(I);
dcoh = adummyvar(S.sig_coh(I));
neural = v(I);
indepvar = {'coh',dcoh,'neural',neural,'session',adummyvar(S.session(I))};
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);

beta(idx.neural)
stats.p(idx.neural)

%% stats for Natalie, with low-coherences only

I = S.choice==0 & RT>min_RT & abs(S.sig_coh)<=0.065;

depvar = RT(I);
dcoh = adummyvar(S.sig_coh(I));
neural = v(I);
indepvar = {'coh',dcoh,'neural',neural};
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);

beta(idx.neural)
stats.p(idx.neural)


