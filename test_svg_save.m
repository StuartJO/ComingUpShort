FitType = 'maxKS';

%mkdir(['./figures/Plot1/',FitType])

mdldata = load('Hansen_networks.mat');

A_dist = mdldata.A_dist;

A = mdldata.adj{1};

d = triu2vec(mdldata.A_dist,1);

avec = triu2vec(A,1);

dist_thr{2} = d<30;
dist_thr{3} = d>=30 & d<=90;
dist_thr{4}= d>90;
dist_thr{1}= d>=0;

d_short = d<30;
d_mid = d>=30 & d<=90;
d_long = d>90;

mdls = 0:9;

NMdls = length(mdls);

MdlFit = zeros(NMdls,100);
DegCorr = zeros(NMdls,100);

MdlDeg = zeros(NMdls,200);

cmap = [0.5 0.5 0.5; lines(7); 0.6941 0.3490 0.1569; [251,154,153]./255];
cmap_alpha = make_alpha_rgb(cmap,.5);

    for form = 1
        for law = 1

            if form == 0
                AddMult = 'Mult';
            else
            AddMult = 'Add';
            end
            if law == 0
                LAW = 'powerlaw';
            else
            LAW = 'exponential';
            end
            form_name = ['PhysMdls_',AddMult,'_',LAW];

            FalseDiscoveryRate = zeros(NMdls,100,4);

for mdlIND = 1:NMdls
    mdl = mdls(mdlIND);
if mdl == 9
    Output = load(['GNM_TopoMdl',num2str(2),'_',AddMult,'_',LAW,'.mat']);
else
Output = load(['GNM_Mdl',num2str(mdl),'_',AddMult,'_',LAW,'.mat']);
end

%load('./data/random200_data4topomdl.mat')
% 
    MdlFit(mdlIND,:) = Output.BestFit.(FitType).maxKS;
    bnets = Output.BestFit.(FitType).b;
    DegCorr(mdlIND,:) = Output.BestFit.(FitType).DegCorr;



n = length(A);

for i = 1:length(bnets)
    b = bnets{i};
    B = zeros(n);
    B(b) = 1;
    B = B + B';
    bvec = triu2vec(B,1);
    MdlDeg(mdlIND,:) = MdlDeg(mdlIND,:)+sum(B);

    for j = 1:4
        thr = dist_thr{j};
        athr = avec.*thr;
        bthr = bvec.*thr;        

        r0(mdlIND,i,j) = sum((athr==0).*(bthr==0))./sum(athr==0);
        r1(mdlIND,i,j) = sum(athr.*bthr)./sum(athr);

        FalseDiscoveryRate(mdlIND,i,j) = sum(athr==0&bthr==1)./sum(bthr);
    end

end

end


MdlDeg = MdlDeg./100;
MAP_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};

%%
%figure('Position',[1 1 1382 747])

PlotKSden = 0;
Addlegend = 0;

if form == 1 && law == 1
    switch FitType
    case 'maxKS'
FigLbl = {'A','B','C'};
FigLbl2 = 'D';
MorS='';
        case 'degcorr'
FigLbl = {'J','K','L'};
MorS='S';
FigLbl2 = 'D';
    end
elseif form == 0 && law == 0
FigLbl = {'A','B','C'};
MorS='S';
FigLbl2 = 'A';
    elseif form == 0 && law == 1
FigLbl = {'D','E','F'};
MorS='S';
FigLbl2 = 'B';
        elseif form == 1 && law == 0
FigLbl = {'G','H','I'};
MorS='S';
FigLbl2 = 'C';
end


figure('Position',[-2000 -300 1888 1117])

for j = 1:3
subplot(2,3,j)
    if j == 1
    x = squeeze(r1(:,:,1));
    y = MdlFit;
    xlabel_name = 'Connection recovery ({\itR})';
    ylabel_name = 'max({\itKS})';
    elseif j == 2
    x = DegCorr;
    y = MdlFit;
    xlabel_name = 'Degree correlation';
    ylabel_name = 'max({\itKS})';
    elseif j == 3
    y = DegCorr;
    x = squeeze(r1(:,:,1));
    xlabel_name = 'Connection recovery ({\itR})';
    ylabel_name = 'Degree correlation';
    end

Grp = ones(size(x));
for i = 1:10; Grp(i,:)=i; end

scatter(x(:),y(:),100,Grp(:),'filled','MarkerFaceAlpha',.5)
colormap(cmap_alpha)
clim([.5 10.5])


hold on
%s1 = scatter(mean(x,2),mean(y,2),100,cmap,'filled','MarkerEdgeColor',[0 0 0]);
for i = 1:size(x,1)
scatter(mean(x(i,:)),mean(y(i,:)),100,cmap(i,:),'filled','MarkerEdgeColor',[0 0 0]);
end

ylabel(ylabel_name)
xlabel(xlabel_name)

set(gca,'FontSize',20)
end

subplot(2,3,4:6)

data = cell(NMdls*3,1);
for i = 1:NMdls
data{i} = squeeze(r1(i,:,2));
data{i+NMdls} = squeeze(r1(i,:,3));
data{i+(NMdls*2)} = squeeze(r1(i,:,4));
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

set(gca,'Position',[0.1300    0.0574    0.7750    0.3981])

lgd = legend(ss,MAP_names,'NumColumns',5,'Location','northoutside');
lgd.FontSize=18;

set(gca,'Position',[0.1300    0.0574    0.7750    0.3981])


AddLetters2Plots(gcf, {'A','B','C','D'},'HShift', -.23, 'VShift', -.07,'FontSize',36)



lgd.Position = [0.0452313280555853 0.0182036466809212 0.917563756086334 0.123198850010589];
lgd.Box = 'off'; 



if Addlegend == 1
lgd = legend(s1,MAP_names,'NumColumns',3, 'Location','southoutside');
lgd.FontSize=12;
lgd.Position = [0.0452313280555853 0.0182036466809212 0.917563756086334 0.123198850010589];
lgd.Box = 'off'; 
end

if PlotKSden == 0
    ax = gca;
ylabel(ylabel_name)
xlabel(xlabel_name)
set(gca,'FontSize',20)

ax.Position = [0.2486    0.1875    0.6564    0.7375];

AddLetters2Plots(gcf, {FigLbl{j}},'HShift', -.23, 'VShift', -.07,'FontSize',36)

end

print(['./figures/Plot1/',FitType,'/',MorS,num2str(j),'_',form_name,'.png'],'-dpng','-r300')


Spos = gcf().Position;

%%


%figure('Position',[110 100 2194 514])
figure('Position',[1 100 1384 514])
%subplot(1,3,1:2)

data = cell(NMdls*3,1);
for i = 1:NMdls
data{i} = squeeze(r1(i,:,2));
data{i+NMdls} = squeeze(r1(i,:,3));
data{i+(NMdls*2)} = squeeze(r1(i,:,4));
end
jittercamp = repmat(cmap,3,1);
JitterPlot(data,jittercamp,1)
xticks([(1+NMdls)/2 ((NMdls+1)+(NMdls*2))/2 ((NMdls*2+1)+(NMdls*3))/2])
xticklabels({'Short-range (<30mm)','Mid-range (30-90mm)','Long-range (>90mm)'})
clear ss
for i = 1:NMdls
    ss(i) = scatter(-1,-1,50,cmap(i,:),'filled');
end
xlim([0.5 (NMdls*3)+.5])
hold on
ylimits = ylim;
plot([NMdls+.5 NMdls+.5],[0 1],'k','LineWidth',2)
plot([(NMdls*2)+.5 (NMdls*2)+.5],[0 1],'k','LineWidth',2)

%leg = legend(ss,MAP_names,'Orientation','vertical','NumColumns',1,'FontSize',10,'Location','northeast');
% leg = legend(ss,MAP_names,'Orientation', 'Horizontal','Location','northoutside','FontSize',10,'NumColumns',5);
% 
% set(gca,'Position',[0.1300    0.1100    0.7750    0.7382])
% axPos = [0.1300    0.1100    0.7750    0.7382];
% 
% aXRange = axPos(3)-axPos(1);
% legXrange = aXRange*.75;
% legXStart = (1-legXrange)/2;
% 
% leg.Position = [legXStart,0.90,legXrange,0.0741];

ylabel({'Proportion of empirical','connections captured'})
set(gca,'FontSize',20)
ax = gca;
%AddLetters2Plots(gcf, {FigLbl2},'HShift', -0.1, 'VShift', -0.1,'FontSize',36)
AddLetters2Plots(gcf, {FigLbl2},'HShift',-.1, 'VShift', -.05,'FontSize',36)

% ax = gca;
% ax.LineWidth = 2;

%print(['./figures/DistThr_',form_name,'.png'],'-dpng','-r300')

% 250% image scaling = 2.5 times figure resolution
% print(['./DistThr_',form_name,'_r0.png'],'-dpng','-r0')
% 
% print(['./DistThr_',form_name,'_r50.png'],'-dpng','-r50')
% 
% print(['./DistThr_',form_name,'_r100.png'],'-dpng','-r100')

% Calculation for final screen resolution = round((R/96)*gcf.Position(3)
% gcf.Position(4))
Spos1 = gcf().Position;

ScatterPlotWidth = round((300/96)*Spos(3));

DesiredWidth = ScatterPlotWidth*3;

SaveRes = round(96*(DesiredWidth/Spos1(3)));

print(['./figures/Plot1/',FitType,'/DistThr_',form_name,'.png'],'-dpng',['-r',num2str(SaveRes)])


%%


figure('Position',[1 100 1384 514])

data = cell(NMdls*4,1);
for i = 1:NMdls
data{i} = squeeze(FalseDiscoveryRate(i,:,1));
data{i+NMdls} = squeeze(FalseDiscoveryRate(i,:,2));
data{i+(NMdls*2)} = squeeze(FalseDiscoveryRate(i,:,3));
data{i+(NMdls*3)} = squeeze(FalseDiscoveryRate(i,:,4));
end
jittercamp = repmat(cmap,4,1);
JitterPlot(data,jittercamp,1)
xticks([(1+NMdls)/2 ((NMdls+1)+(NMdls*2))/2 ((NMdls*2+1)+(NMdls*3))/2 ((NMdls*3+1)+(NMdls*4))/2])
xticklabels({'Overall','Short-range (<30mm)','Mid-range (30-90mm)','Long-range (>90mm)'})
clear ss
for i = 1:NMdls
    ss(i) = scatter(-1,-1,50,cmap(i,:),'filled');
end
xlim([0.5 (NMdls*4)+.5])
hold on
ylimits = ylim;
plot([NMdls+.5 NMdls+.5],[0 1],'k','LineWidth',2)
plot([(NMdls*2)+.5 (NMdls*2)+.5],[0 1],'k','LineWidth',2)
plot([(NMdls*3)+.5 (NMdls*3)+.5],[0 1],'k','LineWidth',2)


ylabel({'Connection ','false discovery rate'})
set(gca,'FontSize',18)

Spos1 = gcf().Position;

ScatterPlotWidth = round((300/96)*Spos(3));

DesiredWidth = ScatterPlotWidth*3;

SaveRes = round(96*(DesiredWidth/Spos1(3)));

print(['./figures/Plot1/',FitType,'/FDRConn_',form_name,'.png'],'-dpng',['-r',num2str(SaveRes)])


%%

figure('Position',[1 100 1384 59])
clear ss
for i = 1:NMdls
    hold on
    %ss(i) = scatter(-1,-1,100,cmap(i,:),'filled');
    ss(i) = plot(nan, nan,'Color','none','MarkerSize', 15,'Marker','o','MarkerFaceColor',cmap(i,:));
end
leg = legend(ss,MAP_names,'Orientation','Horizontal','Location','northoutside','FontSize',16,'NumColumns',5,'Box','off');
box off
axis off
%leg.Position=[-0.0056    0.2741    1.0001    0.4906];
leg.Position=[-0.0056    0.0338    1.0001    0.9713];
print('./figures/Plot1/LEGEND.png','-dpng','-r300')



        end
        end

  
