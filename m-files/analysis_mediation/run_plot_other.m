load for_plots_mediation_norm_to_se.mat

%% choose what to plot
signal = 'TinC'; 
mediator = 'ramp';

%% example plot
m = {leverage.mediator};
s = {leverage.signal};
I = ismember(s,signal) & ismember(m,mediator);

L = leverage(I);

colors = cbrewer('qual','Paired',6);
colors = colors([2,4,6],:);

t = L.t';

p = publish_plot(2,1);
p.next();
[~,h(1)] = niceBars2(t,L.RT.unmediated.m,L.RT.unmediated.se,colors(1,:),0.3);
[~,h(2)] = niceBars2(t,L.RT.mediated_self.m,L.RT.mediated_self.se,colors(2,:),0.5);
[~,h(3)] = niceBars2(t,L.RT.mediated_other.m,L.RT.mediated_other.se,colors(3,:),0.5);
ylabel('Leverage on RT (\rho)');

hl = legend(h,'unmediated','mediated self','mediated other');
set(hl,'location','best');

p.next();
[~,h(1)] = niceBars2(t,L.choice.unmediated.m,L.choice.unmediated.se,colors(1,:),0.3);
[~,h(2)] = niceBars2(t,L.choice.mediated_self.m,L.choice.mediated_self.se,colors(2,:),0.5);
[~,h(3)] = niceBars2(t,L.choice.mediated_other.m,L.choice.mediated_other.se,colors(3,:),0.5);

ylabel('Leverage on choice (\beta)');
xlabel('Time [s]');

%% all vs all
uni_s = unique({leverage.signal});
labels = replace(uni_s,'w_','');

m = {leverage.mediator};
s = {leverage.signal};

X = nan(length(uni_s));
for i=1:length(uni_s)
    for j=1:length(uni_s)
        
        if i==j
            mediation = 1;
        else
        
            I = ismember(s,uni_s(i)) & ismember(m,uni_s(j));
            L = leverage(I);
            t = single(L.t);
            tind1 = findclose(t,0.4);
%             tind2 = findclose(t,0.4);

            mediation = (L.choice.unmediated.m(tind1) - L.choice.mediated_other.m(tind1)) / ...
                (L.choice.unmediated.m(tind1) - L.choice.mediated_self.m(tind1));
            
%             mediation = sum(L.choice.unmediated.m(tind) - L.choice.mediated_other.m(tind))./ ...
%                 sum(L.choice.unmediated.m(tind) - L.choice.mediated_self.m(tind));
        end
        X(i,j) = mediation;
    end
end

p = publish_plot(1,1);
set(gcf,'Position',[610  406  484  295]);
ind = 3:length(uni_s);
imagesc(X(ind,:));
set(gca,'xtick',1:length(uni_s),'xticklabel',labels,'tickdir','out');
set(gca,'ytick',1:length(uni_s),'yticklabel',labels(ind));
colores = cbrewer('seq','YlOrRd',100);
colormap(colores);
colorbar
xlabel('Mediator');
ylabel('Signal');

p.format();
%p.append_to_pdf('fig_degree_mediation',1,1);


%% NSt addition: choice & RT


uni_s = unique({leverage.signal});
labels = replace(uni_s,'w_','');

m = {leverage.mediator};
s = {leverage.signal};

X = nan(length(uni_s));
for i = 1 : length(uni_s)
    for j = 1 : length(uni_s)
        
        if i==j
            mediation = 1;
            mediationRT = 1;
        else
        
            I = ismember(s,uni_s(i)) & ismember(m,uni_s(j));
            L = leverage(I);
            t = single(L.t);
            tind1 = findclose(t,0.4);
%             tind2 = findclose(t,0.4);

            mediation = (L.choice.unmediated.m(tind1) - L.choice.mediated_other.m(tind1)) / ...
                (L.choice.unmediated.m(tind1) - L.choice.mediated_self.m(tind1));
            
            mediationRT = (L.RT.unmediated.m(tind1) - L.RT.mediated_other.m(tind1)) / ...
                (L.RT.unmediated.m(tind1) - L.RT.mediated_self.m(tind1));
            
%             mediation = sum(L.choice.unmediated.m(tind) - L.choice.mediated_other.m(tind))./ ...
%                 sum(L.choice.unmediated.m(tind) - L.choice.mediated_self.m(tind));
        end
        X(i,j) = mediation;
        Xrt(i,j) = mediationRT;
    end
end

indY = [6 5 4 3 7];
indX = [7 3 4 5 6 1 2];

p = publish_plot(1,1);
set(gcf,'Position',[100  100  530  700]);
% ind = 3:length(uni_s);
subplot(2, 1, 1); hold on
imagesc(X(indY, indX));
set(gca,'xtick',1:length(uni_s),'xticklabel',labels(indX),'tickdir','out');
set(gca,'ytick',1:length(uni_s),'yticklabel',labels(indY));
colores = cbrewer('seq','YlOrRd',100);
colormap(colores);
colorbar
axis tight
xlabel('Mediator');
ylabel('Signal');

subplot(2, 1, 2); hold on
imagesc(Xrt(indY, indX));
set(gca,'xtick',1:length(uni_s),'xticklabel',labels(indX),'tickdir','out');
set(gca,'ytick',1:length(uni_s),'yticklabel',labels(indY));
colores = cbrewer('seq','YlOrRd',100);
colormap(colores);
colorbar
axis tight
xlabel('Mediator');
ylabel('Signal');


%% NSt addition: choice & RT pure mediation


uni_s = unique({leverage.signal});
labels = replace(uni_s,'w_','');

m = {leverage.mediator};
s = {leverage.signal};

X = nan(length(uni_s));
for i = 1 : length(uni_s)
    for j = 1 : length(uni_s)
        
        
        if i == j
            I = find(ismember(s,uni_s(i)));
            L = leverage(I(1));
            t = single(L.t);
            tind1 = findclose(t,0.4);
            
            mediation = 1 - (L.choice.mediated_self.m(tind1) / L.choice.unmediated.m(tind1));
            
            mediationRT = 1 - (L.RT.mediated_self.m(tind1) / L.RT.unmediated.m(tind1));
        else
            I = ismember(s,uni_s(i)) & ismember(m,uni_s(j));
            L = leverage(I);
            t = single(L.t);
            tind1 = findclose(t,0.4);
            %             tind2 = findclose(t,0.4);
            
            mediation = 1 - (L.choice.mediated_other.m(tind1) / L.choice.unmediated.m(tind1));
            
            mediationRT = 1 - (L.RT.mediated_other.m(tind1) / L.RT.unmediated.m(tind1));
            
        end
        
        X(i,j) = mediation;
        Xrt(i,j) = mediationRT;
    end
end

indY = [6 5 4 3 7];
indX = [7 3 4 5 6];

p = publish_plot(1,1);
set(gcf,'Position',[100  100  430  700]);
% ind = 3:length(uni_s);
subplot(2, 1, 1); hold on
imagesc(X(indY, indX), [0 1]);
set(gca,'xtick',1:length(uni_s),'xticklabel',labels(indX),'tickdir','out');
set(gca,'ytick',1:length(uni_s),'yticklabel',labels(indY));
colores = cbrewer('seq','YlOrRd',100);
colormap(colores);
colorbar
axis tight
xlabel('Mediator');
ylabel('Signal');

subplot(2, 1, 2); hold on
imagesc(Xrt(indY, indX), [0 1]);
set(gca,'xtick',1:length(uni_s),'xticklabel',labels(indX),'tickdir','out');
set(gca,'ytick',1:length(uni_s),'yticklabel',labels(indY));
colores = cbrewer('seq','YlOrRd',100);
colormap(colores);
colorbar
axis tight
xlabel('Mediator');
ylabel('Signal');




%% single sessions

