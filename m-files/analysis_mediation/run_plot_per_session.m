load for_plots_mediation_norm_to_se.mat

%%

J = ismember({leverage.signal},'MinC');
J = find(J,1);
L1 = leverage(J).per_session.self_mediated;

J = ismember({leverage.signal},'MinI');
J = find(J,1);
L2 = leverage(J).per_session.self_mediated;


p = publish_plot(2,1);
set(gcf,'Position',[556  205  288  607]);

V = {'rho_choice','rho_RT'};

for k = 1:2
    var = V{k};
    p.next();

    n = length(L1);
    for i=1:n
        hold all
        plot(L1(i).tt,L1(i).(var),'b');
    end
    hold all
    h(1) = plot(L1(i).tt,nanmean(cat(2,L1.(var)),2),'b','LineWidth',2);

    n = length(L2);
    for i=1:n
        hold all
        plot(L2(i).tt,L2(i).(var),'r');
    end
    h(2) = plot(L2(i).tt,nanmean(cat(2,L2.(var)),2),'r','LineWidth',2);
    
    symmetric_y(gca);
end

legend(h,'MinC','MinI');

p.current_ax(1);
ylabel('Leverage on choice (\rho)');

p.current_ax(2);
ylabel('Leverage on RT (\beta)');
xlabel('Time from motion onset [s]');

set(p.h_ax,'FontSize',14);
% p.format();

