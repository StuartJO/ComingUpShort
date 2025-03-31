%mkdir(['./figures/Plot1/',FitType])

FitType = 'maxKS';

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

MAP_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};

PlotLables = {'A','B','C'};
ITER = 1;
    for form = 0:1
        for law = 0:1

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


if form == 1 && law == 1
figure('Position',[132 76 1103 700])

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

scatter(x(:),y(:),50,Grp(:),'filled','MarkerFaceAlpha',.5)
colormap(cmap_alpha)
clim([.5 10.5])

hold on

for i = 1:size(x,1)
scatter(mean(x(i,:)),mean(y(i,:)),50,cmap(i,:),'filled','MarkerEdgeColor',[0 0 0]);
end

ylabel(ylabel_name)
xlabel(xlabel_name)

set(gca,'FontSize',12)
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
set(gca,'FontSize',12)

set(gca,'Position',[0.1300    0.0574    0.7750    0.3981])

lgd = legend(ss,MAP_names,'NumColumns',5,'Location','northoutside');
lgd.FontSize=10;

set(gca,'Position',[0.1300    0.0574    0.7750    0.3566])
lgdPosition = lgd.Position;

lgd.Position(1) = ( 1-lgdPosition(3) )/2;

AddLetters2Plots(gcf, {'A','B','C','D'},'HShift', -.06, 'VShift', -0.04,'FontSize',24)

saveas(gcf,['TEST.svg'])


end


if form ~= 1 || law ~= 1

figure('Position',[1 100 1206 514])

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

leg = legend(ss,MAP_names,'Orientation', 'Horizontal','Location','northoutside','FontSize',12,'NumColumns',5);

set(gca,'Position',[0.1300    0.1100    0.7750    0.7382])
axPos = [0.1300    0.1100    0.7750    0.7382];

aXRange = axPos(3)-axPos(1);
legXrange = aXRange*.75;
legXStart = (1-legXrange)/2;

leg.Position = [legXStart,0.90,legXrange,0.0741];

ylabel({'Proportion of empirical','connections captured'})
set(gca,'FontSize',20)
ax = gca;
ax.Position = [0.1009    0.0903    0.8925    0.7658];
leg.Position = [0.0810    0.8938    0.9142    0.0866];
%AddLetters2Plots(gcf, {FigLbl2},'HShift', -0.1, 'VShift', -0.1,'FontSize',36)
AddLetters2Plots(gcf, {PlotLables{ITER}},'HShift',-.1, 'VShift', -.14,'FontSize',36)

print()
%saveas(gcf,['LEN_R1_',PlotLables{ITER},'_',form_name,'.svg'])

ITER = ITER + 1;
end

        end
    end