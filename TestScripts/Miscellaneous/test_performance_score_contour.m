%% Cleanup
close all;
fclose all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(7);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;



%% 
quantitiesOfInterestLowerBound = [1 1];
quantitiesOfInterestUpperBound = [4 4];

quantityOfInterestLowerLimit = [2 2];
quantityOfInterestUpperLimit = [3 3];

nPoint = 1000;


%% individual deficits
quantityOfInterest1 = linspace(quantitiesOfInterestLowerBound(1),quantitiesOfInterestUpperBound(1),nPoint)';
quantityOfInterest2 = linspace(quantitiesOfInterestLowerBound(2),quantitiesOfInterestUpperBound(2),nPoint)';

deficit1 = design_measure_to_deficit(quantityOfInterest1,quantityOfInterestLowerLimit(1),quantityOfInterestUpperLimit(1));

figure;
hold all;
plot(quantityOfInterest1,deficit1,'DisplayName','Deficit');
plot([quantityOfInterestLowerLimit(1) quantityOfInterestLowerLimit(1)],ylim,'k:','Linewidth',0.7,'HandleVisibility','off');
plot([quantityOfInterestUpperLimit(1) quantityOfInterestUpperLimit(1)],ylim,'k:','Linewidth',0.7,'DisplayName','Limits');
grid minor;
plot([quantitiesOfInterestLowerBound(1) quantitiesOfInterestUpperBound(1)],[0 0],'k--','Linewidth',0.7,'HandleVisibility','off');
xlim([quantitiesOfInterestLowerBound(1) quantitiesOfInterestUpperBound(1)]);
xlabel('$$y_j$$','interpreter','latex','FontSize',16);
ylabel('Deficit $$\varphi_j$$','interpreter','latex','FontSize',16);
set(gca,'XTick',[]);
legend;
save_print_figure(gcf,[saveFolder,'PerformanceDeficitOutputSpace'],'Size',[figureSize(1) figureSize(2)/2],'PrintFormat',{'png','pdf'});


%% performance score
[quantityOfInterest1Grid,quantityOfInterest2Grid] = meshgrid(quantityOfInterest1,quantityOfInterest2);
deficitGrid = design_measure_to_deficit([quantityOfInterest1Grid(:),quantityOfInterest2Grid(:)],quantityOfInterestLowerLimit,quantityOfInterestUpperLimit);

[~,performanceScore] = design_deficit_to_label_score(deficitGrid);
performanceScore = reshape(performanceScore, size(quantityOfInterest1Grid));

figure;
hold all;
plot([quantityOfInterestLowerLimit(1) quantityOfInterestLowerLimit(1)],[quantitiesOfInterestLowerBound(2) quantitiesOfInterestUpperBound(2)],'k:','Linewidth',0.7,'HandleVisibility','off');
plot([quantityOfInterestUpperLimit(1) quantityOfInterestUpperLimit(1)],[quantitiesOfInterestLowerBound(2) quantitiesOfInterestUpperBound(2)],'k:','Linewidth',0.7,'HandleVisibility','off');
plot([quantitiesOfInterestLowerBound(1) quantitiesOfInterestUpperBound(1)],[quantityOfInterestLowerLimit(2) quantityOfInterestLowerLimit(2)],'k:','Linewidth',0.7,'HandleVisibility','off');
plot([quantitiesOfInterestLowerBound(1) quantitiesOfInterestUpperBound(1)],[quantityOfInterestUpperLimit(2) quantityOfInterestUpperLimit(2)],'k:','Linewidth',0.7,'DisplayName','Limits');
contour(quantityOfInterest1Grid,quantityOfInterest2Grid,performanceScore,-0.5:0.1:1.5,'DisplayName','Performance Score');
clim([-0.5 1.5]);
colormap(color_palette_tol(9:31,'smooth-rainbow'));
colorbar;
grid minor;
xlabel('$$y_1$$','interpreter','latex','FontSize',16);
ylabel('$$y_2$$','interpreter','latex','FontSize',16);
set(gca,'XTick',[],'YTick',[]);
axis([quantitiesOfInterestLowerBound(1) quantitiesOfInterestUpperBound(1) quantitiesOfInterestLowerBound(2) quantitiesOfInterestUpperBound(2)]);
% title('Performance Score','FontSize',16);
legend;
save_print_figure(gcf,[saveFolder,'PerformanceScoreOutputSpace'],'Size',figureSize,'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

