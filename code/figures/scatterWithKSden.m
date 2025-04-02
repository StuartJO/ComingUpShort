function [PlotMain,PlotTop,PlotSide,Cbar]= scatterWithKSden(x,y,varargin)

p = inputParser;

addRequired(p,'x',@(x) ismatrix(x))
addRequired(p,'y',@(x) ismatrix(x))

addOptional(p,'z',ones(size(x)),@(x) ismatrix(x))

addParameter(p,'xlabel','',@(x)isstring(x))
addParameter(p,'ylabel','',@(x)isstring(x))
addParameter(p,'zlabel','',@(x)isstring(x))

addParameter(p,'Yscale','Linear',@(x)isstring(x)|ischar(x))

addParameter(p,'PlotColorbar','on',@(x)isstring(x)|ischar(x))

addParameter(p,'grouping',ones(size(x)),@(x)ismatrix(x))
addParameter(p,'grouping_colors',[],@(x)ismatrix(x)&size(x,2)==3)

addParameter(p,'colormap',parula,@(x)ismatrix(x)&size(x,2)==3)

addParameter(p,'Annot','',@(x)isstring(x)|ischar(x))

parse(p,x,y,varargin{:});

Xlbl = p.Results.xlabel;
Ylbl = p.Results.ylabel;
Zlbl = p.Results.ylabel;
Yscale = p.Results.Yscale;
Annot = p.Results.Annot;

z = p.Results.z;

grouping = p.Results.grouping;
grouping_colors = p.Results.grouping_colors;
cmap = p.Results.colormap;

plotColorbar = strcmp(p.Results.PlotColorbar,'on');


if strcmp(Annot,'')
PlotAnnot = false;
else
PlotAnnot = true;
end

figure('Position',[100   100   706   694])
mainPlot = subplot(4,4,[5 6 7 9 10 11 13 14 15]);

scatter(x,y,50,z,'filled','MarkerFaceAlpha',.5) 
colormap(cmap)

x_max = max(x,[],'All');
x_min = min(x,[],'All');
ksdenvals_x = linspace(x_min,x_max,100);

if strcmp(Yscale,'log')
y = log(y);
end
y_max = max(y,[],'All');
y_min = min(y,[],'All');

ksdenvals_y = linspace(y_min,y_max,100);

ylabel(Ylbl)
xlabel(Xlbl)

set(gca,'FontSize',16)
xlim([x_min x_max])
xlimits = xlim;
ylimits = ylim;

set(gca, 'LineWidth',2)

set(gca,'YScale',Yscale)
%set(gca, 'YMinorTick','off')

NewYPos=mainPlot.Position(2);
NewYHeight=mainPlot.Position(4);
CurrentYHeight = NewYHeight;
if length(unique(z))~=1 && plotColorbar
%c = colorbar('Position',[0.7871    0.7230    0.0342    0.2644]);
c = colorbar('SouthOutside');
c.Label.String = Zlbl;

NewYPos = mainPlot.Position(2);
NewYHeight = CurrentYHeight-.1;

c.Position = [0.078753534743055,0.081268011527377,0.858923519081308,0.032853021090358];

% NewYPos = 0.13;
% CurrentYPos = mainPlot.Position(2);
% mainPlot.Position(2) = NewYPos;
% YPosDiff = NewYPos-CurrentYPos;
else
    c = [];

end

Grps = unique(grouping);
NGrps = length(Grps);
if isempty(grouping_colors)
    if NGrps <= 7
        grouping_colors = lines(NGrps);
    else
        grouping_colors = turbo(NGrps);
    end
end

smoothKSden = 0;

%%

s_top = subplot(4,4,1:3);

for s = 1:NGrps
    hold on
    xs = x(grouping==Grps(s));
    [U,~,ic] = unique(xs);
    U_counts = accumarray(ic,1);
    range = max(xs)-min(xs);
    bins = 10;
    if smoothKSden
        [denOUT,xout] =  ksdensity(U,ksdenvals_x,'Weights',U_counts); 
    else
        [denOUT,xout] =  ksdensity(xs,ksdenvals_x); 
    end
    plot(xout,denOUT./sum(denOUT),'Color',grouping_colors(s,:),'LineWidth',2)
    %plot(xout,denOUT/bins*range,'Color',grouping_colors(s,:),'LineWidth',2)
end

xlim(xlimits)
ylabel('Proportion')
set(gca,'FontSize',16)
pos = get(gca,'Position');
set(gca,'Position',[pos(1) 0.8350 pos(3) pos(4)])
set(gca,'XAxisLocation','top')
set(gca, 'LineWidth',2)

%s_top.Position(2) = 0.74+YPosDiff;

%%
s_side = subplot(4,4,[8 12 16]);

for s = 1:NGrps
    hold on
    ys = y(grouping==Grps(s));
    [U,~,ic] = unique(ys);
    U_counts = accumarray(ic,1);
    range = max(ys)-min(ys);
    bins = 10;
    if smoothKSden
        [denOUT,yout] =  ksdensity(U,ksdenvals_y,'Weights',U_counts); 
    else
        [denOUT,yout] =  ksdensity(ys,ksdenvals_y); 
    end
    
    if strcmp(Yscale,'log')
    yout = exp(yout);
    end

    plot(denOUT./sum(denOUT),yout,'Color',grouping_colors(s,:),'LineWidth',2)
    %plot(denOUT/bins*range,yout,'Color',grouping_colors(s,:),'LineWidth',2)
end

ylim(ylimits)

xlabel('Proportion')
set(gca,'FontSize',16)
set(gca,'YAxisLocation','right')

set(gca,'YScale',Yscale)

%set(gca, 'YMinorTick','off')

set(gca, 'LineWidth',2)

%%

mainPlot.Position(2) = NewYPos;
mainPlot.Position(4) = NewYHeight;

s_top.Position(2) = NewYPos+NewYHeight+.05;

s_side.Position(2) = NewYPos;
s_side.Position(4) = NewYHeight;

XposStart = mainPlot.Position(1);

mainPlot.Position(1) = XposStart+.01;
s_top.Position(1) = XposStart+.01;

if PlotAnnot
    annot = annotation(gcf, 'textbox',...
        [0,  .92, 0.0 0.0],...
        'String',Annot,...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 36, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";
   
    if plotColorbar
    s_top.Position(2) = NewYPos+NewYHeight+.025;
    end
end

if nargout>0
    PlotMain = mainPlot;
    PlotSide = s_side;
    PlotTop = s_top;
    Cbar = c;
end