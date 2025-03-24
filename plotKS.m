function KS = plotKS(subdata,MdlData,Type,Colors,KScolor)

% for j = 1:length(SubNets)
%     
%     SubNet = double(SubNets{j}>0);
%     
%     subdata{j,1} = sum(SubNet);
%     subdata{j,3} = betweenness_bin(SubNet);
%     subdata{j,2} = clustering_coef_bu(SubNet);
%     subdata{j,4} = Dists(triu(SubNet,1) > 0);    
%     
% end

% figure('Position',[506 536 1539 419])

x_labels = {'Degree, \itk','Clustering, \itc','Betweenness, \itb','Connection distance, \ite'};
y_labels = {'eCDF(\itk)','eCDF(\itc)','eCDF(\itb)','eCDF(\ite)'};

for i = Type

    if iscell(subdata)

    vals = unique([subdata{i} MdlData{i}]);
    for j = 1:length(vals)
        S(j) = mean(vals(j)>=subdata{i});
        S2(j) = mean(vals(j)>=MdlData{i});
    end

    else

   vals = unique([subdata MdlData]);
    for j = 1:length(vals)
        S(j) = mean(vals(j)>=subdata);
        S2(j) = mean(vals(j)>=MdlData);
    end
    


    end


    [KS,I] = max(abs(S2-S));
    Y1 = S2(I);
    Y2 = S(I);
    plot(vals,S2,'Color',Colors,'LineWidth',2)
    hold on
    plot(vals,S,'Color',[0 0 0],'LineWidth',2)
    xlimits = xlim;
    cla
    fill([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[Y1 Y2 Y2 Y1],KScolor,'EdgeColor','none','FaceAlpha',.5)
    hold on
    plot(vals,S2,'Color',Colors,'LineWidth',2)
    hold on
    plot(vals,S,'Color',[0 0 0],'LineWidth',2)
        xlabel(x_labels{i})
    ylabel(y_labels{i})
    box on
end