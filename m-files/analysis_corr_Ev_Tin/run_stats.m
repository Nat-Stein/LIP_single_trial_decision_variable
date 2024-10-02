addpath('../generic/');
addpath('../analysis_reclassif/');

%%

[t, TIN, TOUT, EC, EI, CHOICE, RT, COH, SESSION] = get_and_prep_data();

%%

X = {EC-EI};
Y = {TIN};
XLABEL = {'Time, M_{in}^{contra}-M_{in}^{ipsi} [s]'};
YLABEL = {'Time, T_{in} [s]'};


for k = 1:length(X)
    p = publish_plot(1,1);
    set(gcf,'Position',[531  572  380  320]);
    % colores = cbrewer('div','RdBu',100);
    % colores = cbrewer('div','RdYlBu',100);
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
%             x = EC(:,tind(i)) - EI(:,tind(i));
%             y = TIN(:,tind(j));
            x = X{k}(:,tind(i));
            y = Y{k}(:,tind(j));
            K = ~isnan(x) & ~isnan(y) & I;
%             rho(i,j) = corr2(x(K),y(K));
            [rho(i,j),pval(i,j)] = corr(x(K),y(K));
        end
    end
    
    lim = max(abs(rho(:)));
    lim = [-lim,lim];
    imagesc(t(tind),t(tind),rho',lim);

    colormap(colores);
    colorbar
    axis xy
    h = refline(1,0);
    set(h,'color','k','LineStyle','--');
    xlabel(XLABEL{k});
    ylabel(YLABEL{k});
    p.format('FontSize',14);
    
    %shuffled
    nshuffles = 200;
    rho_shuffle = nan(nt,nt,nshuffles);
    for m=1:nshuffles
        disp(num2str(m));
        xx = X{k};
        yy = Y{k};
        % shuffle rows
        idx = randperm(size(xx,1));
        xx = xx(idx,:);
        for i=1:nt
            for j=1:nt
                x = xx(:,tind(i));
                y = yy(:,tind(j));
                K = ~isnan(x) & ~isnan(y) & I;
                rho_shuffle(i,j,m) = corr(x(K),y(K));
            end
        end
    end


    % make the stats comparison
    tt = t(tind);
    [t1,t2] = meshgrid(tt,tt); %t1: time in ev int dimension; t2: time in momentary ev dimension
    Ja = t1(:)>t2(:) & t1(:)>=0.2 & t2(:)>=0.1;
    
    Jb = t1(:)<t2(:) & t1(:)>=0.2 & t2(:)>=0.1; 
    media = nanmean(rho(Ja))-nanmean(rho(Jb));
    
    % calc on shuffles
    media_shuffle = nan(nshuffles,1);
    for m=1:nshuffles
        r = rho_shuffle(:,:,m);
        media_shuffle(m) = nanmean(r(Ja))-nanmean(r(Jb));
    end
    % fit Gaussian to the shuffles
    pd = fitdist(media_shuffle,'Normal');
    % prob of obtaining a sample as extreme as the one observed
    p_as_extreme(1) = 1-pd.cdf(media);


    % same, but for the other comparison
    
    Jb = t1(:)<=0.2 & t2(:)<=0.1; 
    media = nanmean(rho(Ja))-nanmean(rho(Jb));

    % calc on shuffles
    media_shuffle = nan(nshuffles,1);
    for m=1:nshuffles
        r = rho_shuffle(:,:,m);
        media_shuffle(m) = nanmean(r(Ja))-nanmean(r(Jb));
    end
    % fit Gaussian to the shuffles
    pd = fitdist(media_shuffle,'Normal');
    % prob of obtaining a sample as extreme as the one observed
    p_as_extreme(2) = 1-pd.cdf(media);


%     p.append_to_pdf('fig_corr_after_grouping_across_sessions',k==1,do_save_flag);

end

%%

