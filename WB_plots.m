FitType = 'degcorr';

savdir = ['./figures/WB_plots/',FitType];

mkdir(savdir)

mdldata = load('Hansen_networks_WB.mat');

A_dist = mdldata.A_dist;

A = mdldata.adj{1};

d = triu2vec(A_dist,1);

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

MdlFit = zeros(NMdls,100,1);
DegCorr = zeros(NMdls,100,1);
EdgeOverlap = zeros(NMdls,100,4);
EdgeProp = zeros(NMdls,100,4);
Edge_pdf = zeros(100,length(0:160),NMdls);

MdlEdges = zeros(NMdls,length(d));

f = zeros(NMdls,length(d));
Eta = zeros(NMdls,1);
Gam = zeros(NMdls,1);
Alpha = zeros(NMdls,1);

Nnodes = length(A_dist);

InterHemi = [zeros(Nnodes/2) ones(Nnodes/2); ones(Nnodes/2) zeros(Nnodes/2)]; 

InterHemiVec = triu2vec(InterHemi);

MdlDeg = zeros(NMdls,Nnodes);

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
for mdlIND = 1:NMdls
    mdl = mdls(mdlIND);
if mdl == 8
Output = load(['WB_Hansen_Timing_0_RandMdl',num2str(1),'_',AddMult,'_',LAW,'.mat']);
elseif mdl == 9
    Output = load(['WB_Hansen_Timing_0_Matching_',AddMult,'_',LAW,'.mat']);
else
Output = load(['WB_Hansen_Timing_0_Mdl',num2str(mdl),'_',AddMult,'_',LAW,'.mat']);
end

[~,I] = min(Output.maxKS); 

Output.optim_P = Output.P(I,:);

Eta(mdlIND,:) = Output.optim_P(1);
if length(Output.optim_P)==4
Gam(mdlIND,:) = Output.optim_P(2);
Alpha(mdlIND,:) = Output.optim_P(4);
%f(mdlIND,:) = triu2vec(Output.PDMs{2},1);
elseif length(Output.optim_P)==2
Gam(mdlIND,:) = Output.optim_P(2);
Alpha(mdlIND,:) = 1;
%f(mdlIND,:) = triu2vec(Output.PDMs{2},1);
end 
% 
switch FitType
    case 'maxKS'
    MdlFit(mdlIND,:) = Output.optim_maxKS;
    bnets = Output.optim_b;
    DegCorr(mdlIND,:) = Output.optim_DegCorr;
    case 'degcorr'
    MdlFit(mdlIND,:) = Output.bestDegCorr_maxKS;
    bnets = Output.bestDegCorr_b;
    DegCorr(mdlIND,:) = Output.bestDegCorr_DegCorr;
end

n = length(A);

for i = 1:length(bnets)
    b = bnets{i};
    B = zeros(n);
    B(b) = 1;
    B = B + B';
    bvec = triu2vec(B,1);
    MdlDeg(mdlIND,:) = MdlDeg(mdlIND,:)+sum(B);
    MdlEdges(mdlIND,:) = MdlEdges(mdlIND,:)+bvec';
    [Edge_pdf(i,:,mdlIND),xi] = ksdensity(d(bvec==1),0:160);    
    Edge_cdf(i,:,mdlIND) = ksdensity(d(bvec==1),0:160,'Function','cdf');   
for j = 1:4
        thr = dist_thr{j};
        athr = avec.*thr;
        bthr = bvec.*thr;        
        EdgeOverlap(mdlIND,i,j) = sum(athr & bthr)/sum(athr | bthr);   
        %EdgeOverlap(mdlIND,i,j) = sum(athr.*bthr)./sum(athr);
        %EdgeProp(mdlIND,i,j) = sum(bthr)./sum(thr);
        r0(mdlIND,i,j) = sum((athr==0).*(bthr==0))./sum(athr==0);
        r1(mdlIND,i,j,1) = sum(athr.*bthr)./sum(athr);
        RR(mdlIND,i,j) = sqrt(r0(mdlIND,i,j).*r1(mdlIND,i,j));


        PropCaptured(mdlIND,i,j,1) = sum(athr.*bthr)./sum(athr);
        PropCaptured(mdlIND,i,j,2) = sum(athr&bthr& InterHemiVec==1)./sum(athr & InterHemiVec==1);
        PropCaptured(mdlIND,i,j,3) = sum(athr&bthr& InterHemiVec==0)./sum(athr & InterHemiVec==0);


        FalsePositiveRate(mdlIND,i,j,1) = sum(athr==0&bthr==1)./sum(bthr);
        FalsePositiveRate(mdlIND,i,j,2) = sum(athr==0&bthr==1& InterHemiVec==1)./sum(bthr & InterHemiVec==1);
        FalsePositiveRate(mdlIND,i,j,3) = sum(athr==0&bthr==1& InterHemiVec==0)./sum(bthr & InterHemiVec==0);
end

end

end

MdlEdges= MdlEdges./100;
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

for j = 1:3
figure
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
cmap = [0.5 0.5 0.5; lines(7); 0.6941 0.3490 0.1569; [251,154,153]./255];
cmap_alpha = make_alpha_rgb(cmap,.5);

if PlotKSden == 1
    [MainPlot,~,~,Cbar] = scatterWithKSden(x(:),y(:),Grp(:),'colormap',cmap_alpha,'grouping',Grp(:),'grouping_colors',cmap);
    set(gcf, 'currentaxes', MainPlot)
    delete(Cbar)
else
    scatter(x(:),y(:),100,Grp(:),'filled','MarkerFaceAlpha',.5)
    colormap(cmap_alpha)
    clim([.5 10.5])
end

hold on
%s1 = scatter(mean(x,2),mean(y,2),100,cmap,'filled','MarkerEdgeColor',[0 0 0]);
for i = 1:size(x,1)
s1(i) = scatter(mean(x(i,:)),mean(y(i,:)),100,cmap(i,:),'filled','MarkerEdgeColor',[0 0 0]);
end

if Addlegend == 1
lgd = legend(s1,MAP_names,'NumColumns',3);
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

print([savdir,'/',MorS,num2str(j),'_',form_name,'.png'],'-dpng','-r300')

end

Spos = gcf().Position;
%%
FigLabels = {'A','B','C'};
for ConType = 1:3

figure('Position',[1 100 1384 514])

data = cell(NMdls*4,1);
for i = 1:NMdls
data{i} = squeeze(PropCaptured(i,:,1,ConType));
data{i+NMdls} = squeeze(PropCaptured(i,:,2,ConType));
data{i+(NMdls*2)} = squeeze(PropCaptured(i,:,3,ConType));
data{i+(NMdls*3)} = squeeze(PropCaptured(i,:,4,ConType));
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

if ConType == 1
ylabel({'Proportion of empirical','connections captured'})
elseif ConType == 2
ylabel({'Proportion of empirical','inter-hemispheric','connections captured'})
elseif ConType == 3
ylabel({'Proportion of empirical','intra-hemispheric','connections captured'})
end
set(gca,'FontSize',18)
ax = gca;

AddLetters2Plots(gcf, {FigLabels{ConType}},'HShift',-.11, 'VShift', -.05,'FontSize',36)

Spos1 = gcf().Position;

ScatterPlotWidth = round((300/96)*Spos(3));

DesiredWidth = ScatterPlotWidth*3;

SaveRes = round(96*(DesiredWidth/Spos1(3)));

print([savdir,'/DistThr_',form_name,'_ConType',num2str(ConType),'.png'],'-dpng',['-r',num2str(SaveRes)])

end


%%
FigLabels = {'A','B','C'};
for ConType = 1:3

figure('Position',[1 100 1384 514])

data = cell(NMdls*4,1);
for i = 1:NMdls
data{i} = squeeze(FalsePositiveRate(i,:,1,ConType));
data{i+NMdls} = squeeze(FalsePositiveRate(i,:,2,ConType));
data{i+(NMdls*2)} = squeeze(FalsePositiveRate(i,:,3,ConType));
data{i+(NMdls*3)} = squeeze(FalsePositiveRate(i,:,4,ConType));
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

if ConType == 1
ylabel({'Connection false discovery rate'})
elseif ConType == 2
ylabel({'Inter-hemispheric connection','false discovery rate'})
elseif ConType == 3
ylabel({'Intra-hemispheric connection','false discovery rate'})
end
set(gca,'FontSize',18)
ax = gca;

AddLetters2Plots(gcf, {FigLabels{ConType}},'HShift',-.11, 'VShift', -.05,'FontSize',36)

Spos1 = gcf().Position;

ScatterPlotWidth = round((300/96)*Spos(3));

DesiredWidth = ScatterPlotWidth*3;

SaveRes = round(96*(DesiredWidth/Spos1(3)));

print([savdir,'/FP_DistThr_',form_name,'_ConType',num2str(ConType),'.png'],'-dpng',['-r',num2str(SaveRes)])

end

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

print([savdir,'/DistThr_',form_name,'.png'],'-dpng',['-r',num2str(SaveRes)])

%system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" montage \( ./S1_',form_name,'.png ./S2_',form_name,'.png ./S3_',form_name,'.png \) ./DistThr_',form_name,'.png -geometry +1+2 -tile 1x2 ',form_name,'.png'])

%system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert \( ./S1_',form_name,'.png ./S2_',form_name,'.png ./S3_',form_name,'.png +append \) ./DistThr_',form_name,'.png -append ',form_name,'.png'])

% system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert ( ./figures/Plot1/',FitType,'/S1_',form_name,'.png ./figures/Plot1/,',FitType,'/S2_',form_name,'.png ./figures/Plot1/,',FitType,'/S3_',form_name,'.png +append ) ./figures/Plot1/,',FitType,'/DistThr_',form_name,'.png -append ',form_name,'.png']);
% delete([savdir,'/S1_',form_name,'.png'])
% delete([savdir,'/S2_',form_name,'.png'])
% delete([savdir,'/S3_',form_name,'.png'])
% delete([savdir,'/DistThr_',form_name,'.png'])

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
print('./figures/WB_Plots/LEGEND.png','-dpng','-r300')

run = 0;
if run == 1

load('fsaverage_surface_data.mat')

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

Deg = sum(A);

parc = Scha7_parcs.lh_scha400;
figure('Position',[150.6000  211.4000  560.0000  392.4000])
axes('Position',[0.0029   -0.0972    0.9936    1.2256])
brainCmap = brewermap(256,'YlOrRd');

p_left = plotSurfaceROIBoundary(surface,parc,Deg,'midpoint',brainCmap,2,[0 max(Deg)]);
camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis image

print('./figures/Deg_Lat_0.png','-dpng')

view([90 0])
print('./figures/Deg_Med_0.png','-dpng')

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" montage ./figures/Deg_Lat_',num2str(0),'.png ./figures/Deg_Med_',num2str(0),'.png -geometry +2+1 -tile 2x1 Deg',num2str(0),'.png'])

delete('./figures/Deg_Lat_0.png')
delete('./figures/Deg_Med_0.png')
% 
% system(['"H:\ImageMagick-7.0.10-Q16-HDRI/magick.exe" montage ./figures/OHBMfigures/Deg_Lat_',num2str(0),'.png ./figures/OHBMfigures/Deg_Med_',num2str(0),'.png -geometry +2+1 -tile 2x1 Deg',num2str(0),'.png'])
% 
for i = 1:NMdls

DegData = MdlDeg(i,:);

FaceVertexCData = makeFaceVertexCData(lh_inflated_verts,lh_faces,parc,DegData,brainCmap,[0 max(DegData)],0);

set(p_left,'FaceVertexCData',FaceVertexCData,'EdgeColor','none','FaceColor','flat','Clipping','off');

view([-90 0])
print(['./figures/Deg_Lat_',num2str(i),'.png'],'-dpng')

view([90 0])
print(['./figures/Deg_Med_',num2str(i),'.png'],'-dpng')

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" montage ./figures/Deg_Lat_',num2str(i),'.png ./figures/Deg_Med_',num2str(i),'.png -geometry +2+1 -tile 2x1 ./figures/Deg_',strrep(MAP_names{i}, ' ', ''),'_',form_name,'.png'])

delete(['./figures/Deg_Lat_',num2str(i),'.png'])
delete(['./figures/Deg_Med_',num2str(i),'.png'])

end

end

run = 1;

if run == 1

figure
clear ksplot
for i = 1:NMdls
ksplot(i) = plot(0:160,mean(squeeze(Edge_pdf(:,:,i))),'Color',cmap(i,:),'LineWidth',4);
hold on
end
hold on
Edge_pdf_emp = ksdensity(d(avec==1),0:160);
ksplot_emp = plot(0:160,Edge_pdf_emp,'Color','k','LineWidth',4);
xlabel('Distance (mm)')
ylabel('Proportion of connections')
set(gca,'FontSize',20)
xlim([0 160])
legend([ksplot_emp ksplot],[{'Empirical'} MAP_names],'Orientation','vertical','NumColumns',1,'FontSize',10,'Location','northeast')
print([savdir,'/PDF_',form_name,'.png'],'-dpng')
clf

clear ksplot
for i = 1:NMdls
ksplot(i) = plot(0:160,mean(squeeze(Edge_cdf(:,:,i))),'Color',cmap(i,:),'LineWidth',4);
hold on
end
hold on
Edge_cdf_emp = ksdensity(d(avec==1),0:160,'Function','cdf');
ksplot_emp = plot(0:160,Edge_cdf_emp,'Color','k','LineWidth',4);
xlabel('Distance (mm)')
ylabel('CDF')
set(gca,'FontSize',20)
xlim([0 160])
legend([ksplot_emp ksplot],[{'Empirical'} MAP_names],'Orientation','vertical','NumColumns',1,'FontSize',10,'Location','northeast')

print([savdir,'/CDF_',form_name,'.png'],'-dpng')


%close all
end
        end
    end
