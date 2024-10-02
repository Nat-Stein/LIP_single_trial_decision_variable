    [t,TIN, TOUT,EC,EI,CHOICE,RT,COH,SESSION] = get_and_prep_data();

%%
e = EC-EI;
tind = t>=0.1 & t<=0.4;
ee = nanmean(e(:,tind),2);

%
p = publish_plot(2,1);
p.next();
curva_media(CHOICE, COH, ee<nanmedian(ee), 2);
hold all
curva_media(CHOICE, COH, ee>nanmedian(ee), 2);
hl = legend('low contra-ipsi Min', 'high contra-ipsi Min');

p.next();
curva_media(RT,COH,ee<nanmedian(ee) & CHOICE==0,2);
hold all
curva_media(RT,COH,ee>nanmedian(ee) & CHOICE==0,2);


set(hl,'location','best');

p.format();

p.append_to_pdf('fig_ev_on_behav',1,1);

%% stats
min_RT = 0.5;

% I = ismember(CHOICE,[0,1]) & RT>min_RT;
% depvar = CHOICE(I);
% indepvar = {'coh',adummyvar(COH(I)),'neural',ee(I)};
% [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
% stats.beta(idx.neural)
% stats.p(idx.neural)
% 
% I = CHOICE==0  & RT>min_RT;
% % depvar = zscore_bygroup(RT(I),[COH(I),SESSION(I)]);
% depvar = RT(I);
% indepvar = {'coh',adummyvar(COH(I)),'neural',ee(I)};
% [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
% stats.beta(idx.neural)
% stats.p(idx.neural)

% alternative
I = CHOICE==0  & RT>min_RT & ~isnan(ee);
depvar = zscore_bygroup(RT(I),[COH(I),SESSION(I)]);
% depvar = RT(I);
indepvar = {'neural',ee(I),'ones',ones(sum(I),1)};
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
stats.beta(idx.neural)
stats.p(idx.neural)

% correlation
[rho,pval] = corr(depvar, ee(I));
rho
% pval
tval = (rho * sqrt(sum(I)-2))/(sqrt(1-rho^2)); %t-statistic
pval = tcdf(tval,sum(I)-2)
% P = 2*P; % if two-sided

%% plot the correlation
% I = CHOICE==0 & abs(COH)<0.1;
% depvar = zscore_bygroup(RT(I),[COH(I),SESSION(I)]);

p = publish_plot(1,1);
plot(ee(I),depvar,'.');
lsline
[rho,pval] = corr(ee(I),depvar,'rows','pairwise');
xlabel('Neural activity (MinC-MinI)');
ylabel('RT residuals (z-scored)');
p.text_draw(1, ['p = ',num2str(pval)]);
p.format('MarkerSize',4);
symmetric_x(gca);
symmetric_y(gca);
xlim([-5,5]);
ylim([-5,5]);
axis square
p.append_to_pdf('fig_corr_RT_res',1,1);




%%
% addpath('../analysis_mediation/');
% sProj = EC-EI;
% tind = findclose(t,0.55);
% sProj_mediator = sProj(:,tind);
% 
% % zRT = zscore_bygroup(RT,SESSION);
% 
% minRT = 0.67;
% maxRT = 2;
% start_t = 0;
% end_t = 0.5;
% do_plot = 1;
% flags.norm_to_se = 1;
% I = abs(COH)>0.1;
% [~,out] = corr_with_RT_choice(CHOICE(I)==0, RT(I), COH(I), t, ...
%                         sProj(I,:), sProj_mediator(I), minRT,maxRT, start_t, end_t, do_plot,flags);