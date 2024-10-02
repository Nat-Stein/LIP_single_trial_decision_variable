% load ../data/sigSL_allSessions_Ariel.mat
load ../data/sigSL_allSessions_Ariel_230522.mat

S = sigSL_Ariel;
S.MinI_minus_C = S.MinI - S.MinC;

session = S.session;
times = S.t;
dt = times(2) - times(1);
choice = S.choice;
RT = S.rt/1000;
coh = S.sig_coh;


% str =  {'MinC',
%     'MinI',
%     'PC1',
%     'TinC',
%     'whatD',
%     'whenC',
%     'ramp'};

%%
minRT = 0.67;
maxRT = 2;
nsessions = length(nanunique(session));
N = [];
for j=1:nsessions

    I1 = session==j & abs(coh)<0.1 & RT>minRT & RT<maxRT; % only low coh
    I2 = session==j & abs(coh)<0.1 & RT>minRT; % only low coh
    N = [N; sum(I1), sum(I2)];

end

