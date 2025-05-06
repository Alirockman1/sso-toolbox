%TEST_BATCH_COMPONENT_HOLLOW_SPHERE Batch performance test for hollow sphere
%   TEST_BATCH_COMPONENT_HOLLOW_SPHERE performs component solution space 
%   optimizations in batch, saving the performance data. Factors can be
%	attributed according to the corresponding XLSX file.
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

%% cleanup
close all;
fclose all;
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


%% base 
systemFunction = @distance_to_center;
systemParameter = [0,0,0];
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

% performance requirements
performanceLowerLimit = -inf;
performanceUpperLimit = 5;
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);

% design space 
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
initialPoint = [0,0,0];
componentIndex = {[1,3],[2]};


%% run batch analysis
[solutionSpace,optimizationData,algoData,batchOptions] = batch_sso_stochastic_analysis(...
    'BatchTestHollowSphere.xlsx',...
    designEvaluator,...
    initialPoint,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    componentIndex,...
    'FixedOptions',{'UseAdaptiveGrowthRate',false});


%% plot data
plot_batch_sso_stochastic_analysis_component_metrics(batchOptions,algoData,'SaveFolder',saveFolder);


%% stop transcript 
save([saveFolder,'Data.mat']);
diary off;

