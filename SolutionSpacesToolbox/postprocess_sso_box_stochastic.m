function algorithmData = postprocess_sso_box_stochastic(optimizationData)
%POSTPROCESS_SSO_BOX_STOCHASTIC Postprocess problem and iteration data for box
%   POSTPROCESS_SSO_BOX_STOCHASTIC extracts the most important information from 
%   the data outputs of the SSO stochastic method for boxes and packages it
%   as individual arrays for each iteration step when applicable. This allows 
%   for easy plotting and verification of said data afterwards.
%
%   ALGORITHMDATA = POSTPROCESS_SSO_BOX_STOCHASTIC(PROBLEMDATA,ITERATIONDATA) 
%   receives the data outputs of the SSO stochastic method for boxes in  
%   PROBLEMDATA and ITERATIONDATA and returns the structure with information 
%   arrays in ALGORITHMDATA. PROBLEMDATA should contain all information that is 
%   fixed / the same for all iterations, and ITERATIONDATA should be a structure
%   array with the information that changes each iteration.
%
%   Inputs:
%       - PROBLEMDATA : structure
%       - ITERATIONDATA : (nIter, 1) structure
%
%   Outputs:
%       - ALGORITHMDATA : structure
%           -- IndexExplorationStart : integer
%           -- IndexExplorationEnd : integer
%           -- IndexConsolidationStart : integer
%           -- IndexConsolidationEnd : integer
%           -- IsUsingRequirementSpaces : logical
%           -- GrowthRate : (nIter,1) double
%           -- NumberEvaluatedSamples : (nIter,1) integer
%           -- NumberGoodDesigns : (nIter,1) integer
%           -- NumberPhysicallyFeasibleDesigns : (nIter,1) integer
%           -- NumberAcceptableAndUsefulDesigns : (nIter,1) integer
%           -- NumberAcceptableDesigns : (nIter,1) integer
%           -- NumberUsefulDesigns : (nIter,1) integer
%           -- AlgorithmPhase : (nIter,1) integer
%           -- AlgorithmPhaseIterationNumber : (nIter,1) integer
%           -- MeasureBeforeTrim : (nIter,1) double
%           -- MeasureBeforeTrimNormalized : (nIter,1) double
%           -- MeasureAfterTrim : (nIter,1) double
%           -- MeasureAfterTrimNormalized : (nIter,1) double
%           -- DesignVariableIntervalSizeBeforeTrim : (nIter,nDesignVariable) 
%           double
%           -- DesignVariableIntervalSizeBeforeTrimNormalized : (nIter,
%           nDesignVariable) double
%           -- DesignVariableIntervalSizeAfterTrim : (nIter,nDesignVariable) 
%           double
%           -- DesignVariableIntervalSizeAfterTrimNormalized : (nIter,
%           nDesignVariable) double
%           -- TotalFunctionEvaluations : (nIter,1) integer
%           -- SamplePurity : (nIter,1) double
%           -- RatioMeasureChangeBeforeTrim : (nIter,1) double
%           -- RatioMeasureChangeAfterTrim : (nIter,1) double
%           -- RatioDesignVariableIntervalChangeBeforeTrim : (nIter,
%           nDesignVariable) double
%           -- RatioDesignVariableIntervalChangeAfterTrim : (nIter,
%           nDesignVariable) double
%
%   See also sso_box_stochastic, plot_sso_box_stochastic_metrics, cumsum.
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

    designSpaceBox = [optimizationData.DesignSpaceLowerBound;optimizationData.DesignSpaceUpperBound];
    dsMeasure = optimizationData.Options.MeasureFunction(...
        designSpaceBox,...
        [],...
        optimizationData.Options.MeasureOptions{:});
    designSpaceIntervalSize = optimizationData.DesignSpaceUpperBound - optimizationData.DesignSpaceLowerBound;
    nDesignVariable = size(designSpaceIntervalSize,2);
    
    % algorithm data
    iExplorationStart = 1;
    iConsolidationStart = find([optimizationData.IterationData.Phase]==2,1,'first');
    iExplorationEnd = iConsolidationStart-1;
    iConsolidationEnd = length(optimizationData.IterationData);
    
    % requirement spaces plots only necessary if there are useless designs
    if(all([optimizationData.IterationData.IsUseful],'all'))
        flagReqSpaces = false;
    else
        flagReqSpaces = true;
    end
    
    % create column arrays
    growthRate = [optimizationData.IterationData.GrowthRate]';
    phase = [optimizationData.IterationData.Phase]';
    
    % get measures
    nIter = iConsolidationEnd - iExplorationStart + 1;
    [measureBeforeTrim,...
        measureAfterTrim,...
        nSample,...
        nGood,...
        nPhysicallyFeasible,...
        nAccUse,...
        nAcc,...
        nUse] = deal(nan(nIter,1));
    [designVariableIntervalSizeBeforeTrim,...
        designVariableIntervalSizeAfterTrim] = deal(nan(nIter,nDesignVariable));
    for i=iExplorationStart:iConsolidationEnd
        % measures
        measureBeforeTrim(i) = optimizationData.Options.MeasureFunction(...
            optimizationData.IterationData(i).CandidateBoxBeforeTrim,[],optimizationData.Options.MeasureOptions{:});
        measureAfterTrim(i) = optimizationData.Options.MeasureFunction(...
            optimizationData.IterationData(i).CandidateBoxAfterTrim,[],optimizationData.Options.MeasureOptions{:});

        % number of (labeled) samples
        nSample(i) = size(optimizationData.IterationData(i).EvaluatedDesignSamples,1);
        nGood(i) = sum(optimizationData.IterationData(i).IsGoodPerformance);
        nPhysicallyFeasible(i) = sum(optimizationData.IterationData(i).IsPhysicallyFeasible);
        nAccUse(i) = sum(optimizationData.IterationData(i).IsAcceptable & optimizationData.IterationData(i).IsUseful);
        nAcc(i) = sum(optimizationData.IterationData(i).IsAcceptable);
        nUse(i) = sum(optimizationData.IterationData(i).IsUseful);

        % box interval sizes
        designVariableIntervalSizeBeforeTrim(i,:) = ...
            optimizationData.IterationData(i).CandidateBoxBeforeTrim(2,:) - optimizationData.IterationData(i).CandidateBoxBeforeTrim(1,:);
        designVariableIntervalSizeAfterTrim(i,:) = ...
            optimizationData.IterationData(i).CandidateBoxAfterTrim(2,:) - optimizationData.IterationData(i).CandidateBoxAfterTrim(1,:);
    end

    % iteration number per phase
    phaseIterNumber = [1:iExplorationEnd , ...
        (iConsolidationStart:iConsolidationEnd)-iConsolidationStart+1]';

    % further information
    measureBeforeTrimNormalized = measureBeforeTrim./dsMeasure;
    measureAfterTrimNormalized = measureAfterTrim./dsMeasure;
    designVariableIntervalSizeBeforeTrimNormalized = ...
        designVariableIntervalSizeBeforeTrim./designSpaceIntervalSize;
    designVariableIntervalSizeAfterTrimNormalized = ...
        designVariableIntervalSizeAfterTrim./designSpaceIntervalSize;
    totalEvaluations = cumsum(nSample);
    purity = nAcc./nSample;
    ratioMeasureChangeBeforeTrim = measureBeforeTrim./[0;measureBeforeTrim(1:end-1)];
    ratioMeasureChangeAfterTrim = measureAfterTrim./[0;measureAfterTrim(1:end-1)];
    ratioDesignVariableIntervalChangeBeforeTrim = designVariableIntervalSizeBeforeTrim./...
        [zeros(1,nDesignVariable);designVariableIntervalSizeBeforeTrim(1:end-1,:)];
    ratioDesignVariableIntervalChangeAfterTrim = designVariableIntervalSizeAfterTrim./...
        [zeros(1,nDesignVariable);designVariableIntervalSizeAfterTrim(1:end-1,:)];

    % wrap
    algorithmData = struct(...
    	'IndexExplorationStart',iExplorationStart,...
    	'IndexExplorationEnd',iExplorationEnd,...
    	'IndexConsolidationStart',iConsolidationStart,...
    	'IndexConsolidationEnd',iConsolidationEnd,...
    	'IsUsingRequirementSpaces',flagReqSpaces,...
    	'GrowthRate',growthRate,...
    	'NumberEvaluatedSamples',nSample,...
        'NumberGoodDesigns',nGood,...
        'NumberPhysicallyFeasibleDesigns',nPhysicallyFeasible,...
    	'NumberAcceptableAndUsefulDesigns',nAccUse,...
    	'NumberAcceptableDesigns',nAcc,...
    	'NumberUsefulDesigns',nUse,...
    	'AlgorithmPhase',phase,...
    	'AlgorithmPhaseIterationNumber',phaseIterNumber,...
    	'MeasureBeforeTrim',measureBeforeTrim,...
    	'MeasureBeforeTrimNormalized',measureBeforeTrimNormalized,...
    	'MeasureAfterTrim',measureAfterTrim,...
    	'MeasureAfterTrimNormalized',measureAfterTrimNormalized,...
        'DesignVariableIntervalSizeBeforeTrim',designVariableIntervalSizeBeforeTrim,...
        'DesignVariableIntervalSizeBeforeTrimNormalized',designVariableIntervalSizeBeforeTrimNormalized,...
        'DesignVariableIntervalSizeAfterTrim',designVariableIntervalSizeAfterTrim,...
        'DesignVariableIntervalSizeAfterTrimNormalized',designVariableIntervalSizeAfterTrimNormalized,...
    	'TotalFunctionEvaluations',totalEvaluations,...
    	'SamplePurity',purity,...
    	'RatioMeasureChangeBeforeTrim',ratioMeasureChangeBeforeTrim,...
    	'RatioMeasureChangeAfterTrim',ratioMeasureChangeAfterTrim,...
        'RatioDesignVariableIntervalChangeBeforeTrim',ratioDesignVariableIntervalChangeBeforeTrim,...
        'RatioDesignVariableIntervalChangeAfterTrim',ratioDesignVariableIntervalChangeAfterTrim);
end