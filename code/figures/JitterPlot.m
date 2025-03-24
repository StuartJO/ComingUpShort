function JitterPlot(data,cmap,ColorLine,pos)

jitterOffset = .5;

if nargin < 2 
cmap = lines(length(data));
end


if nargin < 3 
ColorLine=0;
end

if nargin < 4
    pos = 1:length(data);
end


if ~iscell(cmap)
MarkerFaceColor = make_alpha_rgb(cmap,.5);
end

for i = 1:length(data)

jitter = (rand(length(data{i}),1)-.5)*jitterOffset;
x=ones(length(data{i}),1)*pos(i);
%scatter(x+jitter,data{i},'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.1);
if ~iscell(cmap)
scatter(x+jitter,data{i},'filled','MarkerFaceColor',MarkerFaceColor(i,:),'MarkerFaceAlpha',0.75);
else
scatter(x+jitter,data{i},30,cmap{i},'filled','MarkerFaceAlpha',0.75);
end
hold on

plotjitterOffset = jitterOffset/2;
plotjitterOffset = jitterOffset;
if ColorLine==1 && ~iscell(cmap)
plot([pos(i)-plotjitterOffset pos(i)+plotjitterOffset],[nanmean(data{i}) nanmean(data{i})],'LineWidth',2,'Color',cmap(i,:))
else
plot([pos(i)-plotjitterOffset pos(i)+plotjitterOffset],[nanmean(data{i}) nanmean(data{i})],'LineWidth',2,'Color','k')
end
end