
ptiles = [5 10 25 50 100];

TopProb = zeros(150,length(ptiles),10);
ANNOTS = {'A','B','C','D','E','F','G','H','I','J'};
for j = 1:150
dthr = j;
for i = 1:10

    dthrlog = d>=dthr-5 & d<=dthr+5;
    athr = avec(dthrlog);

EXISTprob(i,j) = median(prob(i,dthrlog & avec==1));
NOEXISITprob(i,j) = median(prob(i,dthrlog & avec==0));

    probthr = prob(i,dthrlog);

    for k = 1:length(ptiles)
        ptilethr = prctile(probthr,100-ptiles(k));
        TopProb(j,k,i) = mean(athr(probthr>=ptilethr));
    end

end
end

ptile_leg = cell(1,length(ptiles));
for i = 1:length(ptiles)
ptile_leg{i} = [num2str(ptiles(i)),'%'];
end
ptile_leg{i} = 'Null';
%figure('Position',[1 1 1382 747])
for i = 1:10
    %subplot(2,5,i)
    figure
plot(EXISTprob(i,:),'k','LineWidth',3)
hold on
plot(NOEXISITprob(i,:),'r','LineWidth',3)
set(gca,'Yscale','log')
xlabel('Distance (mm)')
ylabel('Median probability ({\itP_{ij}})')
title(MAP_names{i})
set(gca,'FontSize',16)
set(gca, 'LineWidth',1.5)
annot = annotation(gcf, 'textbox',...
        [0,  .9, 0.0 0.0],...
        'String',ANNOTS{i},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 28, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";
print(['ProbMedianMdl_',ANNOTS{i},'.png'],'-dpng')
end

sky_cmap = [sky(length(ptiles)-1); 0 0 0];
figure('Position',[1 1 1382 747])
for i = 1:10
    %subplot(2,5,i)
    figure
    for j = 1:length(ptiles)
plot(squeeze(TopProb(:,j,i)),'Color',sky_cmap(j,:),'LineWidth',2)
hold on
    end
legend(ptile_leg,'FontSize',14)

xlabel('Distance (mm)')
ylabel({'Proportion of empirical','connections in top X%'})
title(MAP_names{i})
set(gca,'FontSize',16)
set(gca, 'LineWidth',1.5)
annot = annotation(gcf, 'textbox',...
        [0,  .9, 0.0 0.0],...
        'String',ANNOTS{i},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 28, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";
print(['ProbTopMdl_',ANNOTS{i},'.png'],'-dpng')
end
close all

AvgProbgr30mm = mean(TopProb(30:150,3,:),'All');