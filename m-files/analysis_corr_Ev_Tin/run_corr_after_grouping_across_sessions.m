% By Ariel Zylberberg

%%
do_save_flag = 1;

[t, TIN, TOUT, EC, EI, CHOICE, RT, COH, SESSION] = get_and_prep_data();

save(fullfile(saveLoc, 'data_for_corr_analysis'), 't', 'TIN', 'TOUT', 'EC', 'EI', 'CHOICE', 'RT', 'SESSION')

%%

X = {EC-EI,EC, EI, EC,TOUT,cumsum(EC-EI,2)};
Y = {TIN, TIN, TIN, EI,TIN, TIN};
XLABEL = {'Time, M_{in}^{contra}-M_{in}^{ipsi} [s]',...
    'Time, M_{in}^{contra} [s]',...
    'Time, M_{in}^{ipsi} [s]',...
    'Time, M_{in}^{contra} [s]',...
    'Time, T_{out} [s]',...
    'Time, cumsum(M_{in}^{contra}-M_{in}^{ipsi}) [s]',...
    };
YLABEL = {'Time, T_{in} [s]',...
    'Time, T_{in} [s]',...
    'Time, T_{in} [s]',...
    'Time, M_{in}^{ipsi} [s]',...
    'Time, T_{in} [s]',...
    'Time, T_{in} [s]'};

xylims = [-0.01 0.51]; 

for k = 1 % : length(X)
    figure('Position', [531  572  380  320]); hold on
    color_scale = 3;
    switch color_scale
        case 1
            colores = cbrewer('div','RdBu',100);
            %         colores = colores(end:-1:1,:);
        case 2
            %         colores = cbrewer('seq','YlOrRd',100);
            colores = cbrewer('seq','YlGnBu',100);
            colores = colores(end:-1:1,:);
        case 3
            colores = cbrewer('div','RdBu',100);
            colores = colores(end:-1:1,:);
    end
    I = RT>0.55;
    tind = find(t>=0.0 & t<=0.5);
    nt = length(tind);
    rho = nan(nt);
    pval = nan(nt);
    for i=1:nt
        for j=1:nt
            x = X{k}(:,tind(i));
            y = Y{k}(:,tind(j));
            K = ~isnan(x) & ~isnan(y) & I;
            [rho(i,j),pval(i,j)] = corr(x(K),y(K));
        end
    end
    switch color_scale
        case 1
            lim = max(abs(rho(:)));
            lim = [-lim,lim];
            imagesc(t(tind),t(tind),rho',lim);
        case 2
            imagesc(t(tind),t(tind),rho');
        case 3
            lim = max(abs(rho(:)));
            lim = [-lim,lim];
            imagesc(t(tind),t(tind),rho',lim);

    end
    colormap(colores);
    colorbar
    xlim(xylims)
    ylim(xylims)
    axis square
    h = refline(1,0);
    set(h,'color','k','LineStyle','--');
    xlim
    xlabel(XLABEL{k});
    ylabel(YLABEL{k});

end


