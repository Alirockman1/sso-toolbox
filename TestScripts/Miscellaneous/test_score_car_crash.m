%
%   Copyright 2025 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%   
%       http://www.apache.org/licenses/LICENSE-2.0
%   
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

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


%% problem setup
systemFunction = @car_crash_2d;
%               m d1c d2c   v0
systemParameter = [2000 0.3 0.3 15.6];
%                       E_rem    a_max order
performanceLowerLimit = [-inf        0  -inf];
performanceUpperLimit = [   0  32*9.81     0];

bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit,...
    'PerformanceNormalizationFactor',[5.3e5 313.92 1.7e6]...
    ...'PerformanceNormalizationFactor',[nan nan nan]...
    );


%% design space
%                            F1     F2   
designSpaceLowerBound = [   1e5    1e5];
designSpaceUpperBound = [   1e6    1e6];


%% sample in design space
designSample = sampling_latin_hypercube([designSpaceLowerBound;designSpaceUpperBound],100000);
[performanceDeficit,physicalFeasibilityDeficit,evaluationOuptut] = designEvaluator.evaluate(designSample);
[~,performanceScore] = design_deficit_to_label_score(performanceDeficit);
[~,iRobust] = min(performanceScore);



%% Plot Score
figure;
goodTriangleX = [405000 627840 183000 405000];
goodTriangleY = [405000 627840 627840 405000];
plot3(goodTriangleX,goodTriangleY,[0 0 0 0],'g');
hold all;
grid minor;
plot3(designSample(:,1),designSample(:,2),performanceScore,'bo','LineWidth',0.01,'MarkerSize',0.1);
plot3([designSample(iRobust,1) designSample(iRobust,1)],...
      [designSample(iRobust,2) designSample(iRobust,2)],...
      [0.01 min(performanceScore)],...
      'Color','m','Marker','.','MarkerSize',24);
xlabel('F_1 [N]');
ylabel('F_2 [N]');
zlabel('Performance Score');
legend({'Border of Good Designs','Performance Score','Most Robust Design'},'location','southeast');
save_print_figure(gcf,[saveFolder,'PerformanceScore'],'Size',figureSize);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

