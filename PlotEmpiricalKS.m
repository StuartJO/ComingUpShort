function PlotEmpiricalKS(WhatX,WhatY,WhatZ,Parc,TractData)

if ischar(TractData)
switch TractData
    case 'iFOD2'
data = load('.\outputs\Schaefer_7net_iFOD2_acpc_lh_str70Thr_fitMetrics.mat');
    case 'FACT'
data = load('.\outputs\Schaefer_7net_FACT2_acpc_lh_str70Thr_fitMetrics.mat');
end
else
data = TractData;
end

SchSize = 100:100:1000;

SchUse = find(ismember(SchSize,Parc));

cmap = turbo(10);

figure('Position',[474   100   706   694])
subplot(4,4,[5 6 7 9 10 11 13 14 15])

clear y x denOUT xout

sub2use = [1:298 300:973];

switch WhatX
    case 'maxKS'
        X = data.maxKS(sub2use,sub2use,SchUse);
        Xlbl = 'max(\itKS\rm)';
    case 'den'
        X = data.DenDiff(sub2use,sub2use,SchUse);
        Xlbl = 'Density difference';
    case 'DegCorr'
        X = data.DegCorr(sub2use,sub2use,SchUse);
        Xlbl = 'Degree correlation';
    case 'edge'
        X = data.EdgeJaccard(sub2use,sub2use,SchUse);
        Xlbl = 'Connection overlap ({\itJ})';
    case 'edge1'
        X = data.EdgeOverlap{2}(sub2use,sub2use,SchUse);
        Xlbl = 'Short-edge overlap ({\itJ})';
    case 'edge2'
        X = data.EdgeOverlap{3}(sub2use,sub2use,SchUse);
        Xlbl = 'Mid-edge overlap ({\itJ})';
    case 'edge3'
        X = data.EdgeOverlap{4}(sub2use,sub2use,SchUse);
        Xlbl = 'Long-edge overlap ({\itJ})';
    case 'rmse'
        X = data.maxRMSE(sub2use,sub2use,SchUse);
        Xlbl = 'max(RMSE)';
    case 'TFdiff'
        X = data.TFdiff(sub2use,sub2use,SchUse);
        Xlbl = 'TFdiff';
    case 'TND'
        X = data.TopoDist(sub2use,sub2use,SchUse);
        Xlbl = 'TND';
    case 'maxRd'
        X = data.ramRd(sub2use,sub2use,SchUse);
        Xlbl = 'maxRd';
end

switch WhatY
    case 'maxKS'
        Y = data.maxKS(sub2use,sub2use,SchUse);
        Ylbl = 'max(\itKS\rm)';
    case 'den'
        Y = data.DenDiff(sub2use,sub2use,SchUse);
        Ylbl = 'Density difference';
    case 'DegCorr'
        Y = data.DegCorr(sub2use,sub2use,SchUse);
        Ylbl = 'Degree correlation';
    case 'edge'
        Y = data.EdgeJaccard(sub2use,sub2use,SchUse);
        Ylbl = 'Connection overlap ({\itJ})';
    case 'edge1'
        Y = data.EdgeOverlap{2}(sub2use,sub2use,SchUse);
        Ylbl = 'Short-edge overlap ({\itJ})';
    case 'edge2'
        Y = data.EdgeOverlap{3}(sub2use,sub2use,SchUse);
        Ylbl = 'Mid-edge overlap ({\itJ})';
    case 'edge3'
        Y = data.EdgeOverlap{4}(sub2use,sub2use,SchUse);
        Ylbl = 'Long-edge overlap ({\itJ})';
        case 'RMSE'
        Y = data.maxRMSE(sub2use,sub2use,SchUse);
        Ylbl = 'max(RMSE)';
    case 'TFdiff'
        Y = data.TFdiff(sub2use,sub2use,SchUse);
        Ylbl = 'TFdiff';
    case 'TND'
        Y = data.TopoDist(sub2use,sub2use,SchUse);
        Ylbl = 'TND';
    case 'maxRd'
        Y = data.ramRd(sub2use,sub2use,SchUse);
        Ylbl = 'maxRd';
end

switch WhatZ
    case 'maxKS'
        Z = data.maxKS(sub2use,sub2use,SchUse);
        Zlbl = 'max(\itKS\rm)';
    case 'den'
        Z = data.DenDiff(sub2use,sub2use,SchUse);
        Zlbl = 'Density difference';
    case 'DegCorr'
        Z = data.DegCorr(sub2use,sub2use,SchUse);
        Zlbl = 'Degree correlation';
    case 'edge'
        Z = data.EdgeOverlap{1}(sub2use,sub2use,SchUse);
        Zlbl = {'Edge overlap','({\itJ})'};
    case 'edge1'
        Z = data.EdgeOverlap{2}(sub2use,sub2use,SchUse);
        Zlbl = {'Short-edge overlap','({\itJ})'};
    case 'edge2'
        Z = data.EdgeOverlap{3}(sub2use,sub2use,SchUse);
        Zlbl = {'Mid-edge overlap','({\itJ})'};
    case 'edge3'
        Z = data.EdgeOverlap{4}(sub2use,sub2use,SchUse);
        Zlbl = {'Long-edge overlap','({\itJ})'};
    case 'parc'
        Zlbl = 'Schaefer parcellation';
    case 'KS'
        Zlbl = 'KS determinant';
        Zlbls = {'Deg.','Clust.','Bet.','Edge Dst.'};
    case 'RMSE'
        Z = data.maxRMSE(sub2use,sub2use,SchUse);
        Zlbl = 'max(RMSE)';
    case 'TFdiff'
        Z = data.TFdiff(sub2use,sub2use,SchUse);
        Zlbl = 'TFdiff';
    case 'TND'
        Z = data.TopoDist(sub2use,sub2use,SchUse);
        Zlbl = 'TND';
    case 'maxRd'
        Z = data.ramRd(sub2use,sub2use,SchUse);
        Zlbl = 'maxRd';
    case 'corr'        
        Zlbl = 'Corr. determinant';
        Zlbls = {'Deg.','Clust.','Bet.','Nodal Dst.'};
end

smoothKSden = 1;

for s = 1:length(SchUse)
y(s,:) = triu2vec(squeeze(Y(:,:,s)),1);    
x(s,:) = triu2vec(squeeze(X(:,:,s)),1);
end
dist_max = max(x,[],'All');
dist_min = min(x,[],'All');

for s = 1:length(SchUse)
switch WhatZ
    case {'den','maxRd','maxKS','edge','TND','TFdiff','RMSE','DegCorr'}
        z = triu2vec(squeeze(Z(:,:,s)),1);
    case 'KS'
        %[~,z] = max(squeeze(Z(:,:,:,SchUse(s))),[],3); 
        [~,ks] = max(squeeze(data.KS(sub2use,sub2use,SchUse(s),:)),[],3);
        z = triu2vec(ks,1);
    case 'corr'
        [~,Z] = max(1-squeeze(data.TopogCorr(sub2use,sub2use,SchUse(s),:)),[],3);
        z = triu2vec(Z,1);
    otherwise
        z = cmap(SchUse(s),:);
end

   hold on
   scatter(x(s,:),y(s,:),50,z,'filled','MarkerFaceAlpha',.5) 
end

   ylabel(Ylbl)
   xlabel(Xlbl)
   
   set(gca,'FontSize',16)

xlim([dist_min dist_max])

xlimits = xlim;
ylimits = ylim;
set(gca, 'YMinorTick','off')

set(gca, 'LineWidth',2)

switch WhatZ
    case {'parc'}
        % clim([.5 length(SchUse)+.5])
        % c.Ticks = 1:length(SchUse);
        % c.TickLabels = SchSize(SchUse);
        % colormap(cmap(SchUse,:))
        lbl = cell(length(SchUse),1);
        for i = 1:length(SchUse)
             hold on
            ss(i) = plot(nan, nan,'Color','none','MarkerSize', 15,'Marker','o','MarkerFaceColor',cmap(SchUse(i),:));
            lbl{i} = num2str(SchSize(SchUse(i)));
        end
        leg = legend(ss,lbl,'Location','northoutside','FontSize',14,'NumColumns',2,'Box','on');
        leg.Position=[0.71671387157683,0.7415,0.262322948985329,0.198847266782601];
        title(leg,{'Schaefer','parcellation'})
        leg.Position=[0.71671387157683,0.76,0.262322948985329,0.198847266782601];
    case {'KS','corr'}
        c = colorbar('Position',[0.7871    0.7230    0.0342    0.2644]);
        c.Label.String = Zlbl;
        clim([.5 4.5])
        c.Ticks = 1:4;
        c.TickLabels = Zlbls;
        colormap(lines(4))
    otherwise
        c = colorbar('Position',[0.7871    0.7230    0.0342    0.2644]);
        c.Label.String = Zlbl;
        clim([min(Z,[],'All') max(Z,[],'All')])
end

  % 
s_top = subplot(4,4,1:3);
ksdenvals_x = linspace(dist_min,dist_max,100);

for s = 1:length(SchUse)
    hold on
    xs = x(s,:);
    [U,~,ic] = unique(xs);
    U_counts = accumarray(ic,1);

    if smoothKSden
        [denOUT,xout] =  ksdensity(U,ksdenvals_x,'Support','positive','Weights',U_counts); 
    else
        [denOUT,xout] =  ksdensity(xs,ksdenvals_x,'Support','positive'); 
    end

plot(xout,denOUT,'Color',cmap(SchUse(s),:),'LineWidth',2)

end

xlim(xlimits)
ylabel('Proportion')
set(gca,'FontSize',16)
pos = get(gca,'Position');
set(gca,'Position',[pos(1) 0.8350 pos(3) pos(4)])

set(gca,'XAxisLocation','top')
  % 
  % ecdf(prob(avec==0))
set(gca, 'LineWidth',2)

s_top.Position(2) = 0.74;

subplot(4,4,[8 12 16])
  
prob_max = max(y,[],'All');
prob_min = min(y,[],'All');
ksdenvals_y = linspace(prob_min,prob_max,100);

for s = 1:length(SchUse)
    hold on
    ys = y(s,:);
    [U,~,ic] = unique(ys);
    U_counts = accumarray(ic,1);

    if smoothKSden
        [denOUT,xout] =  ksdensity(U,ksdenvals_y,'Support','positive','Weights',U_counts); 
    else
        [denOUT,xout] =  ksdensity(y(s,:),ksdenvals_y,'Support','positive'); 
    end

plot(denOUT,xout,'Color',cmap(SchUse(s),:),'LineWidth',2)
end

%set(gca,'YScale','log')
ylim(ylimits)

xlabel('Proportion')
set(gca,'FontSize',16)
set(gca,'YAxisLocation','right')
set(gca, 'YMinorTick','off')

set(gca, 'LineWidth',2)