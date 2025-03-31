function PlotGNMPerformance(INPUT,PLOTLABELS,SAVEDIR)

load(INPUT,'R','EdgeFDR','maxKS','DegCorr','Mdl_names')

if nargin < 2
    PLOTLABELS = {'A','B','C','D','E'};
end

if nargin < 3
    [name, ~] = fileparts(INPUT);
    SAVEDIR=name;
end

mkdir(['./figures/',SAVEDIR])

NMdls = length(Mdl_names);

if NMdls == 10
cmap = [0.5 0.5 0.5; lines(7); 0.6941 0.3490 0.1569; [251,154,153]./255];
elseif NMdls == 12
cmap = [166,206,227;...
31,120,180;...
178,223,138;...
51,160,44;...
251,154,153;...
227,26,28;...
253,191,111;...
255,127,0;...
202,178,214;...
106,61,154;...
255,255,153;...
177,89,40]./255;
end
cmap_alpha = make_alpha_rgb(cmap,.5);

for j = 1:3
figure
    if j == 1
    x = squeeze(R(:,:,1));
    y = maxKS;
    xlabel_name = 'Connection recovery ({\itR})';
    ylabel_name = 'max({\itKS})';
    elseif j == 2
    x = DegCorr;
    y = maxKS;
    xlabel_name = 'Degree correlation';
    ylabel_name = 'max({\itKS})';
    elseif j == 3
    y = maxKS;
    x = squeeze(R(:,:,1));
    xlabel_name = 'Connection recovery ({\itR})';
    ylabel_name = 'Degree correlation';
    end

Grp = ones(size(x));
for i = 1:NMdls; Grp(i,:)=i; end

scatter(x(:),y(:),100,Grp(:),'filled','MarkerFaceAlpha',.5)
colormap(cmap_alpha)
clim([.5 NMdls+.5])

hold on

for i = 1:size(x,1)
scatter(mean(x(i,:)),mean(y(i,:)),100,cmap(i,:),'filled','MarkerEdgeColor',[0 0 0]);
end

ylabel(ylabel_name)
xlabel(xlabel_name)

ax = gca;
ylabel(ylabel_name)
xlabel(xlabel_name)
set(gca,'FontSize',20)

ax.Position = [0.2486    0.1875    0.6564    0.7375];

AddLetters2Plots(gcf, {PLOTLABELS{j}},'HShift', -.23, 'VShift', -.07,'FontSize',36)

print(['./figures/',SAVEDIR,'/Scatter',num2str(j),'.png'],'-dpng','-r300')

end

Spos = gcf().Position;

figure('Position',[1 100 1385 514])

data = cell(NMdls*3,1);
for i = 1:NMdls
data{i} = squeeze(R(i,:,2));
data{i+NMdls} = squeeze(R(i,:,3));
data{i+(NMdls*2)} = squeeze(R(i,:,4));
end
jittercamp = repmat(cmap,3,1);
JitterPlot(data,jittercamp,1)
xticks([(1+NMdls)/2 ((NMdls+1)+(NMdls*2))/2 ((NMdls*2+1)+(NMdls*3))/2])
xticklabels({'Short-range (<30mm)','Mid-range (30-90mm)','Long-range (>90mm)'})
clear ss
for i = 1:NMdls
    ss(i) = scatter(-1,-1,50,cmap(i,:),'filled','MarkerEdgeColor',[0 0 0]);
end
xlim([0.5 (NMdls*3)+.5])
hold on
ylimits = ylim;
plot([NMdls+.5 NMdls+.5],[0 1],'k','LineWidth',2)
plot([(NMdls*2)+.5 (NMdls*2)+.5],[0 1],'k','LineWidth',2)

ylabel({'Proportion of empirical','connections captured'})
set(gca,'FontSize',20)
Spos1 = gcf().Position;
ScatterPlotWidth = round((300/96)*Spos(3));
DesiredWidth = ScatterPlotWidth*3;
SaveRes = round(96*(DesiredWidth/Spos1(3)));

set(gca,'Position',[0.0990    0.1100    0.8815    0.8150])

AddLetters2Plots(gcf, {PLOTLABELS{4}},'HShift',-.09, 'VShift', -.05,'FontSize',36)

print(['./figures/',SAVEDIR,'/DistThr.png'],'-dpng',['-r',num2str(SaveRes)])

%%

figure('Position',[50 100 1384 84])
clear ss
for i = 1:NMdls
    hold on
    %ss(i) = scatter(-1,-1,100,cmap(i,:),'filled');
    ss(i) = plot(nan, nan,'Color','none','MarkerSize', 15,'Marker','o','MarkerFaceColor',cmap(i,:));
end
leg = legend(ss,Mdl_names,'Orientation','Horizontal','Location','northoutside','FontSize',16,'NumColumns',5,'Box','off');
box off
axis off
%leg.Position=[-0.0056    0.2741    1.0001    0.4906];
leg.Position=[-0.0022    0.1870    1.0007    0.6964];

print(['./figures/',SAVEDIR,'/LEGEND.png'],'-dpng',['-r',num2str(SaveRes)])

%%

%system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert \( ./figures/',SAVEDIR,'/Scatter1.png ./figures/',SAVEDIR,'/Scatter2.png ./figures/',SAVEDIR,'/Scatter3.png +append \) ./figures/',SAVEDIR,'/DistThr.png ./figures/',SAVEDIR,'/LEGEND.png -append final.png'])


%system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" montage ./figures/Deg_Lat_',num2str(0),'.png ./figures/Deg_Med_',num2str(0),'.png -geometry +2+1 -tile 2x1 Deg',num2str(0),'.png'])
