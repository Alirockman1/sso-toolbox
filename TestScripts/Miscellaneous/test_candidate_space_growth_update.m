%test_component_hollow_sphere Component solution spaces for a sphere problem 
%   test_component_hollow_sphere computes a component solution spaces with 
%   for a sphere problem.
%
%   Copyright 2024 Eduardo Rodrigues Della Noce
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
RNGstate = rng;
savefolder = save_diary_files(mfilename);
goldenratio = (1+sqrt(5))/2;
figure_size = [goldenratio 1]*8.5;


%% function call
%
systemFunction = @distance_to_center;
systemParameter = [0,0,0];
%        x1  x2  x3
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
componentIndex = {[1,3],[2]};
%
performanceLowerLimit = -inf;
performanceUpperLimit = 5;
%
initialRegion = [-5,-5,-5;5,5,5];
nSample = 50;
trimmingOptions = {...
    'TrimmingMethodFunction',@component_trimming_method_corner_box_removal...
    };


bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);
nComponent = length(componentIndex);
for i=1:nComponent
    candidateSpace(i) = CandidateSpaceCornerBoxRemoval(designSpaceLowerBound(componentIndex{i}),designSpaceUpperBound(componentIndex{i}));
end

% initial sample and trimming
designSampleInitial = sampling_latin_hypercube(initialRegion,nSample);
[performanceDeficit,physicalFeasibilityDeficit] = designEvaluator.evaluate(designSampleInitial);
[isGoodPerformance,performanceScore] = design_deficit_to_label_score(performanceDeficit);
if(~isempty(physicalFeasibilityDeficit))
    isPhysicallyFeasible = design_deficit_to_label_score(physicalFeasibilityDeficit);
else
    isPhysicallyFeasible = true(size(designSampleInitial,1),1);
end
[isAcceptable,isUseful] = design_requirement_spaces_label('Omega2',isGoodPerformance,isPhysicallyFeasible);
trimmingOrder = trimming_order(~isAcceptable,performanceScore);
candidateSpaceInitial = component_trimming_operation(designSampleInitial,isAcceptable&isUseful,...
    trimmingOrder,componentIndex,candidateSpace,trimmingOptions{:});

% growth
for i=1:nComponent
    candidateSpaceGrown(i) = candidateSpaceInitial(i).grow_candidate_space(0.01);
end


%% Plot Comparison
figure;
hold all;
candidateSpaceInitial(1).plot_candidate_space(gcf,'EdgeColor',color_palette_tol('cyan'),'FaceColor','none');
plot(designSampleInitial(isGoodPerformance,1),designSampleInitial(isGoodPerformance,3),'go');
plot(designSampleInitial(~isGoodPerformance,1),designSampleInitial(~isGoodPerformance,3),'rx');

shapeInitial = candidateSpaceInitial(1).DesignSampleDefinition(candidateSpaceInitial(1).IsShapeDefinition,:);
shapeGrown = candidateSpaceGrown(1).DesignSampleDefinition(candidateSpaceGrown(1).IsShapeDefinition,:);

figure;
hold all;
candidateSpaceInitial(1).plot_candidate_space(gcf,'EdgeColor',color_palette_tol('cyan'),'FaceColor','none');
candidateSpaceGrown(1).plot_candidate_space(gcf,'EdgeColor',color_palette_tol('green'),'FaceColor','none');
plot(shapeInitial(:,1),shapeInitial(:,2),'g.');
plot(shapeGrown(:,1),shapeGrown(:,2),'b.');
grid minor;


%% new trim
% new trim
designSampleNew = sampling_latin_hypercube(initialRegion,nSample);
[performanceDeficit,physicalFeasibilityDeficit] = designEvaluator.evaluate(designSampleNew);
[isGoodPerformance,performanceScore] = design_deficit_to_label_score(performanceDeficit);
if(~isempty(physicalFeasibilityDeficit))
    isPhysicallyFeasible = design_deficit_to_label_score(physicalFeasibilityDeficit);
else
    isPhysicallyFeasible = true(size(designSampleNew,1),1);
end
[isAcceptable,isUseful] = design_requirement_spaces_label('Omega2',isGoodPerformance,isPhysicallyFeasible);
trimmingOrder = trimming_order(~isAcceptable,performanceScore);
candidateSpaceTrimmed = component_trimming_operation(designSampleNew,isAcceptable&isUseful,...
    trimmingOrder,componentIndex,candidateSpaceGrown,trimmingOptions{:});

figure;
hold all;
candidateSpaceTrimmed(1).plot_candidate_space(gcf,'EdgeColor',color_palette_tol('cyan'),'FaceColor','none');
plot(designSampleNew(isGoodPerformance,1),designSampleNew(isGoodPerformance,3),'go');
plot(designSampleNew(~isGoodPerformance,1),designSampleNew(~isGoodPerformance,3),'rx');



%% Save and Stop Transcripting
save([savefolder,'Data.mat']);
diary off;

