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
nDimension = 1;
targetPurity = 0.7;
exponent = 1.1;

purity = 0:0.001:1;
growthRateBase = 0.1;
minimumGrowthRate = 0;

growthRate1 = ((1-targetPurity)*purity./((1-purity)*targetPurity)).^(exponent./nDimension).*growthRateBase;
growthRate2 = ((1-targetPurity)./(1-purity)).^(exponent./nDimension).*growthRateBase;

figure;
plot(purity,growthRate1);
grid minor;
hold all;
plot(purity,growthRate2);
plot(purity,purity./targetPurity*growthRateBase);
lim = axis;
plot([0 1],[growthRateBase growthRateBase],'k--');
plot([targetPurity targetPurity],lim(3:4),'k:');
axis([min(purity) max(purity) lim(3) lim(4)]);
xlabel('Purity $$p_{i-1}$$','interpreter','latex');
ylabel('Growth Rate $$g_i$$','interpreter','latex')
legend({'Modified Scheme 1 - $$(1-p_{i-1})/p_{i-1}$$','Modified Scheme 2 - $$1-p_{i-1}$$',...
    'Original Scheme','Previous Growth Rate $$g_{i-1}$$','Target Purity $$p^t$$'},...
    'Location','northwest','interpreter','latex');
lim = axis;
axis([...
    0 ... % x min
    1 ... % x max
    0 ... % y min
    min(lim(4),0.2) % y max
    ]);
% save_print_figure(gcf,[saveFolder,'ModifiedGrowthRate'],'Size',figureSize,'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

