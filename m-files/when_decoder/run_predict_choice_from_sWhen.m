load ../data/sigSL_allSessions_Ariel_230522.mat
S = sigSL_Ariel;
load ../data/sigRL_allSessions_Ariel_230522.mat
R = sigRL_Ariel;

do_save_fig = 0;

%% Natalie will re-do this for the version that's in the paper (06/05/2023)

%%

% av_win = 0.01;
av_win = 0.05;
sm = round(av_win/(S.t(2)-S.t(1)));

% v = 'TinC';
v = 'whenC';

% stim aligned
signal_s = conv2(1,ones(sm,1),S.(v),'same'); % smooth the signal
ss = signal_s - nanmean(signal_s); % remove the mean
A = nansum((ss>0).*(S.choice==0) + (ss<0).*(S.choice==1));
A = A./sum(~isnan(ss));

% resp aligned
signal_r = conv2(1,ones(sm,1),R.(v),'same'); % smooth the signal
ss = signal_r - nanmean(signal_r); % remove the mean
B = nansum((ss>0).*(R.choice==0) + (ss<0).*(R.choice==1));
B = B./sum(~isnan(ss));

%% auc and regression 

% auc
tind = S.t>=0 & S.t<=0.5;
tstim = S.t(tind);
ss = signal_s(:,tind);
for i=1:size(ss,2)
    [~,~,~,auc_stim(i)] = perfcurve(S.choice, ss(:,i),0);
    % regression
    I = ~isnan(ss(:,i));
    indepvar = {'neural',ss(I,i),'ones',ones(sum(I),1)};
    [beta,idx,stats,x,LRT] = f_regression(S.choice(I)==0,[],indepvar);
    yhat = glmval(beta,x,'logit','constant','off');
    regre_stim(i) = mean([(yhat>0.5)==S.choice(I)==0]);

end

%
tind = R.t>=-0.3 & R.t<=0.0;
tresp = R.t(tind);
ss = signal_r(:,tind);
for i=1:size(ss,2)
    [~,~,~,auc_resp(i)] = perfcurve(R.choice, ss(:,i),0);
    % regression
    I = ~isnan(ss(:,i));
    indepvar = {'neural',ss(I,i),'ones',ones(sum(I),1)};
    [beta,idx,stats,x,LRT] = f_regression(R.choice(I)==0,[],indepvar);
    yhat = glmval(beta,x,'logit','constant','off');
    regre_resp(i) = mean([(yhat>0.5)==R.choice(I)==0]);
end

%% save for Natalie
tind = S.t>=0 & S.t<=0.5;
acc_pred.motion_aligned.t = S.t(tind);
acc_pred.motion_aligned.pred = A(tind);
tind = R.t>=-0.3 & R.t<=0;
acc_pred.resp_aligned.t  = R.t(tind);
acc_pred.resp_aligned.pred = B(tind);
save choice_pred_from_Swhen acc_pred


%% plot
p = publish_plot(1,2);
set(gcf,'Position',[302  307  667  374]);
p.shrink(1:2,1,0.8);
p.displace_ax(1:2,0.07,2);
p.next();
plot(S.t,A,'k');
hold on
xlim([0,0.5]);
ylim([0.4,1]);

p.next();
plot(R.t,B,'k');
hold on
xlim([-0.3,0]);
ylim([0.4,1]);



%
p.current_ax(1);
p.shrink(1:2,1,0.8);
p.displace_ax(1:2,0.07,2);
plot(tstim,auc_stim,'r');
plot(tstim,regre_stim,'b');

p.current_ax(2);
plot(tresp,auc_resp,'r');
plot(tresp,regre_resp,'b');


%
p.same_xscale();
same_ylim(p.h_ax);

p.unlabel_center_plots();

p.format();

p.current_ax(1);
xlabel('Time from motion onset [s]');
ylabel('Choice dec. accuracy by Swhen');
legend('sign trick','auc','regression');
p.current_ax(2);
xlabel('Time from resp. [s]');

p.append_to_pdf('fig_choice_from_sWhen',1,do_save_fig);


