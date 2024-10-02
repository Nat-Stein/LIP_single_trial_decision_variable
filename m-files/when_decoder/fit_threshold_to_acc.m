function [p, thres, prop_cross,out] = fit_threshold_to_acc(dv, coh, correct)

v_thres = linspace(1,20,1000);
for ithres=1:length(v_thres)
    thres = v_thres(ithres);
    crossed = any(dv>thres,1)';
    u = unique(abs(coh));

    for i=1:length(u)
        pc(i) = sum(coh==-u(i) & crossed)./sum(abs(coh)==u(i) & crossed);
    end

    pc(u==0) = 0.5;

    [~,xx] = curva_media(correct, abs(coh),[],0);
    
    good(ithres) = RSquared(xx,pc');
    PC(:,ithres) = pc;
    
end

J = any(isnan(PC),1)==0;
I = find(good==max(good(J)),1);
thres = v_thres(I);
prop_cross = mean(any(dv>thres,1));

p = publish_plot(1,1);
plot(u,PC(:,I),'k-');
hold all
[tt,xx,ss] = curva_media(correct, abs(coh),[],0);
errorbar(tt,xx,ss,'Color','k','Marker','o','MarkerFaceColor','k',...
    'MarkerEdgeColor','k','LineStyle','none');
xlabel('Motion strength');
ylabel('Accuracy');


out.acc = xx;
out.acc_model = PC(:,I);
out.ucoh = u;
p.format('FontSize',14,'LineWidthPlot',1);


