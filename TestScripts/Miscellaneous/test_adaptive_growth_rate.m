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
targetPurity = 0.7;
growthExponent = [1,3,9];
colorExponent = color_palette_tol({'blue','red','green'});
purity = 0:0.0001:1;

legendEntry = {};
for i=1:size(growthExponent,2)
    growthRateFactor1(i,:) = ((1-targetPurity)*purity./((1-purity)*targetPurity)).^(1/growthExponent(i));
    growthRateFactor2(i,:) = ((1-targetPurity)./(1-purity)).^(1./growthExponent(i));
    [legendEntry{end+1},legendEntry{end+2}] = deal(...
        sprintf('$$g_e = %d$$ - Update Strategy 1',growthExponent(i)),...
        sprintf('$$g_e = %d$$ - Update Strategy 2',growthExponent(i)));
end

figure;
hold all;
for i=1:size(growthExponent,2)
    plot(purity,growthRateFactor1(i,:),'-','color',colorExponent(i,:),'linewidth',2.5);
    plot(purity,growthRateFactor2(i,:),'-.','color',colorExponent(i,:),'linewidth',2.5);
end
plot(purity,purity./targetPurity,'--','color',color_palette_tol('grey'),'linewidth',2.5);
grid minor;
lim = axis;
plot([0 1],[1 1],'k--','linewidth',1.0,'HandleVisibility','off');
plot([targetPurity targetPurity],lim(3:4),'k:','linewidth',1.5);
axis([min(purity) max(purity) lim(3) lim(4)]);
xlabel('Purity $$a_{i-1}$$','interpreter','latex','FontSize',14);
ylabel('Growth Rate Adaptation Factor $$(g_i/g_{i-1})$$','interpreter','latex','FontSize',14);
legend([legendEntry,{'Original Strategy','Target Purity $$p^t$$'}],...
    'Location','northwest','interpreter','latex','FontSize',10);
lim = axis;
axis([...
    0 ... % x min
    1 ... % x max
    0 ... % y min
    min(lim(4),2) % y max
    ]);
save_print_figure(gcf,[saveFolder,'ModifiedGrowthRate'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

