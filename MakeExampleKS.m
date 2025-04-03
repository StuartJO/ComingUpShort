% This code just makes up some distributions to plot to make examples of
% what CDFs look like and show how maxKS works

% EMP{1} = randi(50,50,1);
% MDL{1} = randi(48,50,1);
% 
% EMP{2} = rescale(randn(50,1));
% MDL{2} = rescale(randn(50,1));
% 
% EMP{3} = rand(50,1);
% MDL{3} = rand(50,1);
% 
% EMP{4} = rescale(exprnd(1,[1000 1]),10,150);
% MDL{4} = rescale(exprnd(6,[1000 1]),10,100);
for i = 1:4
KS(i) = plotKS(EMP,MDL,i,lines(1),[.5 .5 .5]);
xlimits = xlim;
set(gca,'FontSize',24)
xlim(xlimits)
set(gca, 'LineWidth',1.5)
axis square
print(['ExampleKS',num2str(i),'_back.png'],'-dpng','-r300')
close all
plotKS(EMP,MDL,i,lines(1),[1 1 1]);
xlimits = xlim;
set(gca,'FontSize',24)
xlim(xlimits)
set(gca, 'LineWidth',1.5)
print(['./figures/ExampleKS',num2str(i),'.png'],'-dpng','-r300')
close all
end