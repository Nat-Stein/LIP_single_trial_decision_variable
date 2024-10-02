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

str =  {'MinC',
    'MinI',
    'PC1',
    'TinC',
    'whatD',
    'whenC',
    'ramp',
    'MinI_minus_C'};



% smooth in window
for i=1:length(str)
    sm = round(0.05/dt);
    h = ones(sm,1)/sm;
    S.(str{i}) = conv2(1, h, S.(str{i}), 'same');
end

%% time subset
tind = findclose(times, -0.1:0.01:0.75);
for i=1:length(str)
    aux = S.(str{i});
    S.(str{i}) = aux(:,tind);
end
times = times(tind);


%%

flags.norm_to_se = 1;
nsessions =  8;

count = 0;
for iMethod=1:length(str)
    sProj = S.(str{iMethod});
    disp(num2str(iMethod));
    for iMediator=1:length(str)

        if iMediator~=iMethod

            tind = findclose(times, 0.55);
            sProj_mediator = sProj(:,tind);
            sOther_mediator = S.(str{iMediator})(:,tind);

            count = count + 1;

            %% all sessions together

            minRT = 0.67;
            maxRT = 2;
            start_t = 0; % 0
            %     end_t = 0.5;
            end_t = 0.6;

            %% per session and average
            do_plot = 0;
            clear out out_other
            i = 0;
            for j=1:nsessions
                I = session==j & abs(coh)<0.1; % only low coh
                %         I = session==i;
                % i = i+1;
                if ~[all(isnan(to_vec(sProj(I,:)))) || all(isnan(sProj_mediator(I))) || all(isnan(sOther_mediator(I)))]
                    i = i+1;
                    [~,out(i)] = corr_with_RT_choice(choice(I)==0, RT(I), coh(I), times, ...
                        sProj(I,:), sProj_mediator(I), minRT,maxRT, start_t, end_t, do_plot,flags);
                    [~,out_other(i)] = corr_with_RT_choice(choice(I)==0, RT(I), coh(I), times, ...
                        sProj(I,:), sOther_mediator(I), minRT,maxRT, start_t, end_t,do_plot,flags);
                end

            end

            %% plot
            p = publish_plot(2,1);
            p.next();
            [hl(1)] = niceBars2(out(1).tt,nanmean(cat(2,out.rho_RT),2),stderror(cat(2,out.rho_RT),2),'b');
            hold all
            [hl(2)] = niceBars2(out(1).tt,nanmean(cat(2,out.rho_RT_partial),2),stderror(cat(2,out.rho_RT_partial),2),'r');
            [hl(3)] = niceBars2(out(1).tt,nanmean(cat(2,out_other.rho_RT_partial),2),stderror(cat(2,out_other.rho_RT_partial),2),'g');
            ylabel('\rho RT');
            h = legend(hl,'unmediated','mediated by itself','mediated by other');

            p.next();
            niceBars2(out(1).tt,nanmean(cat(2,out.rho_choice),2),stderror(cat(2,out.rho_choice),2),'b');
            hold all
            niceBars2(out(1).tt,nanmean(cat(2,out.rho_choice_partial),2),stderror(cat(2,out.rho_choice_partial),2),'r');
            niceBars2(out(1).tt,nanmean(cat(2,out_other.rho_choice_partial),2),stderror(cat(2,out_other.rho_choice_partial),2),'g');
            xlabel('Time [s]');
            ylabel('\beta choice');

            p.text_draw_fig(['Signal: ',str{iMethod},', Mediator:',str{iMediator}]);

            p.format();
            set(h,'location','best','box','off','fontsize',10);

            if flags.norm_to_se==1
                p.append_to_pdf('fig_mediation_all_methods_norm_to_SE',iMethod==1,1);
            else
                p.append_to_pdf('fig_mediation_all_methods',iMethod==1,1);
            end

            close all

            %% save

            leverage(count).signal = str{iMethod};
            leverage(count).mediator = str{iMediator};

            leverage(count).t = out(1).tt;
            

            x = cat(2, out.rho_RT);
            [leverage(count).RT.unmediated.m , leverage(count).RT.unmediated.se] = averageCorrelation(x);

            x = cat(2, out.rho_RT_partial);
            [leverage(count).RT.mediated_self.m, leverage(count).RT.mediated_self.se] = averageCorrelation(x);

            x = cat(2, out_other.rho_RT_partial);
            [leverage(count).RT.mediated_other.m, leverage(count).RT.mediated_other.se] = averageCorrelation(x);


            leverage(count).choice.unmediated.m = nanmean(cat(2,out.rho_choice),2);
            leverage(count).choice.unmediated.se = stderror(cat(2,out.rho_choice),2);
            leverage(count).choice.mediated_self.m = nanmean(cat(2,out.rho_choice_partial),2);
            leverage(count).choice.mediated_self.se = stderror(cat(2,out.rho_choice_partial),2);
            leverage(count).choice.mediated_other.m = nanmean(cat(2,out_other.rho_choice_partial),2);
            leverage(count).choice.mediated_other.se = stderror(cat(2,out_other.rho_choice_partial),2);
            
            leverage(count).per_session.self_mediated = out;
            leverage(count).per_session.other_mediated = out_other;
            
        end
    end

end

if flags.norm_to_se==1
    save('for_plots_mediation_norm_to_se','leverage');
else
    save('for_plots_mediation','leverage');
end