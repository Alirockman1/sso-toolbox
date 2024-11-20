%test_bounding_box Visualization of the relaxed/strict bounding boxes
%   test_bounding_box allows the visualization of the difference between the
%	variants of the bounding box, namely its strict and relaxed versions.
%
% Copyright 2024 Eduardo Rodrigues Della Noce
% SPDX-License-Identifier: Apache-2.0

% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


%% Cleanup
close all;
fclose all;
clear all;
clc;
more off;
diary off;

%%
designSpaceLowerBound = [-2,-2];
designSpaceUpperBound = [2,2];

designSampleGood = rand(200,2);
designSampleBad = rand(100,2);

designSampleGood = -1 + designSampleGood*2;
designSampleBad = -2 + designSampleBad*4;

designSample = [designSampleGood;designSampleBad];
label = [true(size(designSampleGood,1),1);false(size(designSampleBad,1),1)];

[boundingBoxStrict,boundingBoxRelaxed] = design_bounding_box(designSample,label);

candidateSpace = CandidateSpaceBoundingBox(designSpaceLowerBound,designSpaceUpperBound);
candidateSpace = candidateSpace.define_candidate_space(designSample,label);
shapeDefinition = candidateSpace.DesignSampleDefinition(candidateSpace.IsShapeDefinition,:);

labelBoundary = design_find_boundary_samples(designSample,label);

figure;
plot(designSampleGood(:,1),designSampleGood(:,2),'g.','MarkerSize',10);
hold all;
plot(designSampleBad(:,1),designSampleBad(:,2),'r.','MarkerSize',10);
plot(designSample(labelBoundary,1),designSample(labelBoundary,2),'ko','MarkerSize',10);
plot(shapeDefinition(:,1),shapeDefinition(:,2),'m*','MarkerSize',10);
plot_design_box_2d(gcf,boundingBoxStrict,'EdgeColor','c');
plot_design_box_2d(gcf,boundingBoxRelaxed,'EdgeColor','m');
grid minor;
