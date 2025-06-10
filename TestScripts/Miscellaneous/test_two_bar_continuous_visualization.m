%% Cleanup
fclose all;
close all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(4);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%%
nDivisions = 10;

for i=1:nDivisions
    figure;
    hold all;

    % compute area distribution
    areaDistribution1 = sampling_random([0.06;0.2],i);
    areaDistribution2 = sampling_random([0.06;0.2],i);
    if(i==1)
        tick1 = [0,0.5];
        tick2 = [0.5,1];
        areaDistribution1 = [areaDistribution1,areaDistribution1];
        areaDistribution2 = [areaDistribution2,areaDistribution2];
    else
        tick1 = linspace(0,0.5,i);
        tick2 = linspace(0.5,1,i);
    end

    % smaller divisions
    if(i>2)
        for j=2:i-1
            plot([tick1(j),tick1(j)],[0.2,0],'--','Linewidth',0.5,'Color',color_palette_tol('grey'));
            plot([tick2(j),tick2(j)],[0.2,0],'--','Linewidth',0.5,'Color',color_palette_tol('grey'));
        end
    end

    % truss 
    nodePosition = [0 0; 0.5 0; 1 0];
    nodeElement = [1 2; 2 3];
    plot_truss_deformation(gcf,nodePosition,nodeElement,'TrussPlotOptions',{'Color',color_palette_tol('blue')},'MaximumLinewidth',5.0,'ShowBarNumber',true,'BarNumberOptions',{'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16});

    % walls
    wallX1 = [0 0];
    wallX2 = [1 1];
    wallY = [-0.06 0.06];
    wallOptions = {'Linewidth',8.0,'Color','k'};
    handleWall1 = plot(wallX1,wallY,wallOptions{:});
    handleWall2 = plot(wallX2,wallY,wallOptions{:});

    % force and end-divisions
    plot(0.5,0,'k.','MarkerSize',20);
    quiver(0.5,-0.05,0.2,0,'Color',color_palette_tol('red'),'Linewidth',4.0);
    plot([0.5,0.5],[0.2,-0.05],'k--','Linewidth',2.0);
    plot([0,0],[0.2,-0.05],'k--','Linewidth',2.0);
    plot([1,1],[0.2,-0.05],'k--','Linewidth',2.0);

    % area distribution
    plot(tick1,areaDistribution1,'-o','Color',color_palette_tol('green'),'Linewidth',2.0);
    plot(tick2,areaDistribution2,'-o','Color',color_palette_tol('yellow'),'Linewidth',2.0);

    % text
    if(i==1)
        text(0.25,areaDistribution1(1),'$$A_{(1)}$$','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16,'Color',color_palette_tol('green'),'Interpreter','latex');
        text(0.75,areaDistribution2(1),'$$A_{(2)}$$','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16,'Color',color_palette_tol('yellow'),'Interpreter','latex');
    else
        for j=1:i
            text(tick1(j),areaDistribution1(j),sprintf('$$A_{(1)%d}$$',j),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',16,'Color',color_palette_tol('green'),'Interpreter','latex');
            text(tick2(j),areaDistribution2(j),sprintf('$$A_{(2)%d}$$',j),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16,'Color',color_palette_tol('yellow'),'Interpreter','latex');
        end
    end

    % axis
    axis([0,1,-0.06,0.2]);
    set(gca,'XColor', 'none','YColor','none');

    % save
    save_print_figure(gcf,[saveFolder,sprintf('TwoBarContinuousVisualization-%02d',i)],'Size',figureSize*1.1,'PrintFormat',{'png','pdf'});
end


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;
