mdldata = load('Hansen_networks.mat');

A_dist = mdldata.A_dist;

A = mdldata.adj{1};

C = zeros(160,1);

deg = sum(A);
for i = 1:160
    deg_thr(i,:) = sum(A.*(A_dist>=i));
    C(i) = corr(deg',deg_thr(i,:)');
    %dprop(i,:) = sum(d>=i)/length(d);
end

figure('Position',[125.8000  140.2000  400.8000  538.8000])
plot(C,'Color','k','LineWidth',4)
ylabel('Correlation with empirical degree')
xlabel('Distance threshold')
set(gca,'FontSize',20)
ax = gca;
ax.LineWidth = 2;

%print(['./figures/RevDegThr.png'],'-dpng','-r300')
exportgraphics(gcf,['./figures/RevDegThr.png'],'resolution',300)
% yyaxis right
% plot(dprop.*100,'LineWidth',4)
% ylabel('% connections in threshold')

figure
s = scatterfit(deg',deg_thr(30,:)',100,[227 26 28]./255,' ',1,24);
xlabel('Full empirical degree')
ylabel('Thresholded degree')
xticks(0:10:50)
ax = gca;
ax.LineWidth = 2;
print(['./figures/RevDegThr30Scatter.png'],'-dpng','-r300')

figure
s = scatterfit(deg',deg_thr(60,:)',100,[227 26 28]./255,' ',1,24);
xlabel('Full empirical degree')
ylabel('Thresholded degree')
xticks(0:10:50)
ax = gca;
ax.LineWidth = 2;
print(['./figures/RevDegThr60Scatter.png'],'-dpng','-r300')

figure
s = scatterfit(deg',deg_thr(90,:)',100,[227 26 28]./255,' ',1,24);
xlabel('Full empirical degree')
ylabel('Thresholded degree')
xticks(0:10:50)
ax = gca;
ax.LineWidth = 2;
print(['./figures/RevDegThr90Scatter.png'],'-dpng','-r300')

load('fsaverage_surface_data.mat')

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

Deg = sum(A);

MaxDeg = max(Deg);

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

print('./figures/RevDeg_lat_thr_0.png','-dpng')

view([90 0])
print('./figures/RevDeg_med_thr_0.png','-dpng')

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" montage ./figures/RevDeg_lat_thr_',num2str(0),'.png ./figures/RevDeg_med_thr_',num2str(0),'.png -geometry +2+1 -tile 2x1 ./figures/RevDegthr',num2str(0),'.png'])


for i = [30 60 90]

DegData = deg_thr(i,:);

FaceVertexCData = makeFaceVertexCData(lh_inflated_verts,lh_faces,parc,DegData,brainCmap,[0 MaxDeg],0);

set(p_left,'FaceVertexCData',FaceVertexCData,'EdgeColor','none','FaceColor','flat','Clipping','off');

view([-90 0])
print(['./figures/RevDeg_lat_thr_',num2str(i),'.png'],'-dpng')

view([90 0])
print(['./figures/RevDeg_med_thr_',num2str(i),'.png'],'-dpng')

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" montage ./figures/RevDeg_lat_thr_',num2str(i),'.png ./figures/RevDeg_med_thr_',num2str(i),'.png -geometry +2+1 -tile 2x1 ./figures/RevDegthr',num2str(i),'.png'])

end
