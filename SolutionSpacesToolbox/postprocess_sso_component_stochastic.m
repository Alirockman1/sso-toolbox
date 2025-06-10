function algorithmData = postprocess_sso_component_stochastic(optimizationData)
%POSTPROCESS_SSO_COMPONENT_STOCHASTIC Postprocess for stochastic component SSO
%   POSTPROCESS_SSO_COMPONENT_STOCHASTIC extracts the most important information  
%   from the data outputs of the SSO stochastic method for components and 
%   packages it as individual arrays for each iteration step. This allows for 
%   easy plotting and verification of said data afterwards.
%
%   ALGORITHMDATA = POSTPROCESS_SSO_COMPONENT_STOCHASTIC(PROBLEMDATA, 
%   ITERATIONDATA) receives the data outputs of the SSO stochastic method for
%   components in PROBLEMDATA and ITERATIONDATA and returns the structure with 
%   information arrays in ALGORITHMDATA. PROBLEMDATA should contain all 
%   information that is fixed / the same for all iterations, and ITERATIONDATA 
%   should be a structure array with the information that changes each 
%   iteration.
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
%           -- NumberPaddingSamplesGenerated : (nIter,1) integer
%           -- NumberPaddingSamplesUsed : (nIter,1) integer
%           -- NumberGoodDesigns : (nIter,1) integer
%           -- NumberPhysicallyFeasibleDesigns : (nIter,1) integer
%           -- NumberAcceptableAndUsefulDesigns : (nIter,1) integer
%           -- NumberAcceptableDesigns : (nIter,1) integer
%           -- NumberUsefulDesigns : (nIter,1) integer
%           -- AlgorithmPhase : (nIter,1) integer
%           -- AlgorithmPhaseIterationNumber : (nIter,1) integer
%           -- ComponentMeasureBeforeTrim : (nIter,nComponent) double
%           -- ComponentMeasureBeforeTrimNormalized : (nIter,nComponent) double
%           -- ComponentMeasureAfterTrim : (nIter,nComponent) double
%           -- ComponentMeasureAfterTrimNormalized : (nIter,nComponent) double
%           -- TotalMeasureBeforeTrim : (nIter,1) double
%           -- TotalMeasureBeforeTrimNormalized : (nIter,1) double
%           -- TotalMeasureAfterTrim : (nIter,1) double
%           -- TotalMeasureAfterTrimNormalized : (nIter,1) double
%           -- TotalFunctionEvaluations : (nIter,1) integer
%           -- SamplePurity : (nIter,1) double
%           -- RatioComponentMeasureChangeBeforeTrim : (nIter,nComponent) double
%           -- RatioComponentMeasureChangeAfterTrim : (nIter,nComponent) double
%           -- RatioTotalMeasureChangeBeforeTrim : (nIter,1) double
%           -- RatioTotalMeasureChangeAfterTrim : (nIter,1) double
%           -- RatioPadding : (nIter,1) double
%
%   See also sso_component_stochastic, plot_sso_component_stochastic_metrics.
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

	designSpaceLowerBound = optimizationData.DesignSpaceLowerBound;
    designSpaceUpperBound = optimizationData.DesignSpaceUpperBound;
    designSpaceMeasure = prod(designSpaceUpperBound - designSpaceLowerBound);
    
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
    nPadGenerated = [optimizationData.IterationData.NumberPaddingSamplesGenerated]';
    phase = [optimizationData.IterationData.Phase]';
    timeElapsedAdaptGrowthRate = [optimizationData.IterationData.TimeElapsedAdaptGrowthRate]';
    timeElapsedGrow = [optimizationData.IterationData.TimeElapsedGrow]';
    timeElapsedGenerate = [optimizationData.IterationData.TimeElapsedGenerate]';
    timeElapsedEvaluate = [optimizationData.IterationData.TimeElapsedEvaluate]';
    timeElapsedLabel = [optimizationData.IterationData.TimeElapsedLabel]';
    timeElapsedCount = [optimizationData.IterationData.TimeElapsedCount]';
    timeElapsedShape = [optimizationData.IterationData.TimeElapsedShape]';
    timeElapsedPrepare = [optimizationData.IterationData.TimeElapsedPrepare]';
    timeElapsedTrimmingOrder = [optimizationData.IterationData.TimeElapsedTrimmingOrder]';
    timeElapsedTrim = [optimizationData.IterationData.TimeElapsedTrim]';
    timeElapsedLeanness = [optimizationData.IterationData.TimeElapsedLeanness]';
    timeElapsedMeasure = [optimizationData.IterationData.TimeElapsedMeasure]';
    timeElapsedConvergence = [optimizationData.IterationData.TimeElapsedConvergence]';
    timeElapsedIteration = [optimizationData.IterationData.TimeElapsedIteration]';
    
    % get measures
    componentIndex = optimizationData.ComponentIndex;
    nComponent = size(optimizationData.ComponentIndex,2);
    componentMeasuresBeforeTrim = nan(iConsolidationEnd,nComponent);
    componentMeasuresAfterTrim = nan(iConsolidationEnd,nComponent);
    componentMeasuresBeforeTrimNormalized = nan(iConsolidationEnd,nComponent);
    componentMeasuresAfterTrimNormalized = nan(iConsolidationEnd,nComponent);
    for i=1:iConsolidationEnd
        for j=1:nComponent
            if(isempty(optimizationData.IterationData(i).CandidateSpacesBeforeTrim(j)) || isempty(optimizationData.IterationData(i).CandidateSpacesBeforeTrim(j).DesignSampleDefinition))
                if(isempty(optimizationData.IterationData(i).CandidateSpacesAfterTrim))
                    componentMeasuresBeforeTrim(i,j) = 0;
                else
                    samplingBox = optimizationData.IterationData(i).SamplingBoxBeforeTrim;
                    sbU = samplingBox(2,componentIndex{j});
                    sbL = samplingBox(1,componentIndex{j});
                    componentMeasuresBeforeTrim(i,j) = prod(sbU - sbL);
                end
            else
                componentMeasuresBeforeTrim(i,j) = optimizationData.IterationData(i).ComponentMeasureBeforeTrim(j);
            end

            if(isempty(optimizationData.IterationData(i).CandidateSpacesAfterTrim))
                componentMeasuresAfterTrim(i,j) = componentMeasuresBeforeTrim(i,j);
            else
                componentMeasuresAfterTrim(i,j) = optimizationData.IterationData(i).ComponentMeasureAfterTrim(j);
            end

            componentMeasuresBeforeTrimNormalized(i,j) = componentMeasuresBeforeTrim(i,j)...
            	./prod(designSpaceUpperBound(componentIndex{j}) - designSpaceLowerBound(componentIndex{j}));
            componentMeasuresAfterTrimNormalized(i,j) = componentMeasuresAfterTrim(i,j)...
            	./prod(designSpaceUpperBound(componentIndex{j}) - designSpaceLowerBound(componentIndex{j}));
        end
    end
    measureBeforeTrim = prod(componentMeasuresBeforeTrim,2);
    measureAfterTrim = prod(componentMeasuresAfterTrim,2);
    measureBeforeTrimNormalized = measureBeforeTrim./designSpaceMeasure;
    measureAfterTrimNormalized = measureAfterTrim./designSpaceMeasure;

    % 
    nIter = iConsolidationEnd - iExplorationStart + 1;
    [nSample, nPadUsed, nGood, nPhysicallyFeasible, nAccUse, nAcc, nUse] = ... 
        deal(nan(nIter,1));
    for i=iExplorationStart:iConsolidationEnd
        nSample(i) = size(optimizationData.IterationData(i).EvaluatedDesignSamples,1);
        nPadUsed(i) = size(optimizationData.IterationData(i).PaddingSamplesUsed,1);
        nGood(i) = sum(optimizationData.IterationData(i).IsGoodPerformance);
        nPhysicallyFeasible(i) = sum(optimizationData.IterationData(i).IsPhysicallyFeasible);
        nAccUse(i) = sum(optimizationData.IterationData(i).IsAcceptable & optimizationData.IterationData(i).IsUseful);
        nAcc(i) = sum(optimizationData.IterationData(i).IsAcceptable);
        nUse(i) = sum(optimizationData.IterationData(i).IsUseful);
    end
    phaseIterNumber = [1:iExplorationEnd , ...
        (iConsolidationStart:iConsolidationEnd)-iConsolidationStart+1]';

    % further information
    totalEvaluations = cumsum(nSample);
    purity = nAcc./nSample;
    ratioComponentMeasureChangeBeforeTrim = componentMeasuresBeforeTrim./...
        [zeros(1,nComponent);componentMeasuresBeforeTrim(1:end-1,:)];
    ratioComponentMeasureChangeAfterTrim = componentMeasuresAfterTrim./...
        [zeros(1,nComponent);componentMeasuresAfterTrim(1:end-1,:)];
    ratioTotalMeasureChangeBeforeTrim = measureBeforeTrim./[0;measureBeforeTrim(1:end-1)];
    ratioTotalMeasureChangeAfterTrim = measureAfterTrim./[0;measureAfterTrim(1:end-1)];
    ratioPaddingSamples = nPadUsed./(nPadUsed+nSample);
    totalTimeElapsed = cumsum(timeElapsedIteration);

    % wrap
    algorithmData = struct(...
    	'IndexExplorationStart',iExplorationStart,...
    	'IndexExplorationEnd',iExplorationEnd,...
    	'IndexConsolidationStart',iConsolidationStart,...
    	'IndexConsolidationEnd',iConsolidationEnd,...
    	'IsUsingRequirementSpaces',flagReqSpaces,...
    	'GrowthRate',growthRate,...
    	'NumberEvaluatedSamples',nSample,...
    	'NumberPaddingSamplesGenerated',nPadGenerated,...
    	'NumberPaddingSamplesUsed',nPadUsed,...
        'NumberGoodDesigns',nGood,...
        'NumberPhysicallyFeasibleDesigns',nPhysicallyFeasible,...
    	'NumberAcceptableAndUsefulDesigns',nAccUse,...
    	'NumberAcceptableDesigns',nAcc,...
    	'NumberUsefulDesigns',nUse,...
    	'AlgorithmPhase',phase,...
    	'AlgorithmPhaseIterationNumber',phaseIterNumber,...
    	'ComponentMeasureBeforeTrim',componentMeasuresBeforeTrim,...
    	'ComponentMeasureBeforeTrimNormalized',componentMeasuresBeforeTrimNormalized,...
    	'ComponentMeasureAfterTrim',componentMeasuresAfterTrim,...
    	'ComponentMeasureAfterTrimNormalized',componentMeasuresAfterTrimNormalized,...
    	'TotalMeasureBeforeTrim',measureBeforeTrim,...
    	'TotalMeasureBeforeTrimNormalized',measureBeforeTrimNormalized,...
    	'TotalMeasureAfterTrim',measureAfterTrim,...
    	'TotalMeasureAfterTrimNormalized',measureAfterTrimNormalized,...
    	'TotalFunctionEvaluations',totalEvaluations,...
    	'SamplePurity',purity,...
        'RatioComponentMeasureChangeBeforeTrim',ratioComponentMeasureChangeBeforeTrim,...
        'RatioComponentMeasureChangeAfterTrim',ratioComponentMeasureChangeAfterTrim,...
    	'RatioTotalMeasureChangeBeforeTrim',ratioTotalMeasureChangeBeforeTrim,...
    	'RatioTotalMeasureChangeAfterTrim',ratioTotalMeasureChangeAfterTrim,...
    	'RatioPadding',ratioPaddingSamples,...
        'TotalTimeElapsed',totalTimeElapsed,...
        'TimeElapsedAdaptGrowthRate',timeElapsedAdaptGrowthRate,...
        'TimeElapsedGrow',timeElapsedGrow,...
        'TimeElapsedGenerate',timeElapsedGenerate,...
        'TimeElapsedEvaluate',timeElapsedEvaluate,...
        'TimeElapsedLabel',timeElapsedLabel,...
        'TimeElapsedCount',timeElapsedCount,...
        'TimeElapsedShape',timeElapsedShape,...
        'TimeElapsedPrepare',timeElapsedPrepare,...
        'TimeElapsedTrimmingOrder',timeElapsedTrimmingOrder,...
        'TimeElapsedTrim',timeElapsedTrim,...
        'TimeElapsedLeanness',timeElapsedLeanness,...
        'TimeElapsedMeasure',timeElapsedMeasure,...
        'TimeElapsedConvergence',timeElapsedConvergence,...
        'TimeElapsedIteration',timeElapsedIteration);
end