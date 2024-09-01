function [designSample,paddingSample] = candidate_space_sampling_individual_feasible(candidateSpace,componentIndex,nSample,varargin)
%CANDIDATE_SPACE_SAMPLING_INDIVIDUAL_FEASIBLE Sampling inside candidate spaces
%   CANDIDATE_SPACE_SAMPLING_INDIVIDUAL_FEASIBLE produces design sample points 
%   that are inside all of the candidates spaces; it does so by sequentially 
%   generating samples that are inside each candidate space individually, and 
%   then combining the results at the end. While not enough samples have been 
%   generated for each component separately, the process repeats for those 
%   components.
%
%   DESIGNSAMPLE = CANDIDATE_SPACE_SAMPLING_INDIVIDUAL_FEASIBLE(CANDIDATESPACE,
%   COMPONENTINDEX,NSAMPLE) generates NSAMPLE design sample points that are  
%   inside the candidate spaces CANDIDATESPACE, with components COMPONENTINDEX, 
%   returning said sample in DESIGNSAMPLE.
%
%   DESIGNSAMPLE = CANDIDATE_SPACE_SAMPLING_INDIVIDUAL_FEASIBLE(...NAME,VALUE,
%   ...) allows for the specification of additional options. These are:
%       - 'SamplingMethodFunction' : base sampling method to be used. Default:
%       @sampling_latin_hypercube.
%       - 'SamplingMethodOptions' : extra options for the base sampling method.
%       Default is empty.
%
%   [DESIGNSAMPLE,PADDINGSAMPLE] = CANDIDATE_SPACE_SAMPLING_INDIVIDUAL_FEASIBLE
%   (...) additionally returns extra samples generated PADDINGSAMPLE, which are
%   not inside at least one of the candidate spaces.
%
%   Input:
%       - CANDIDATESPACE : (1,nCandidateSpace) CandidateSpaceBase
%       - COMPONENT : (1,nComponent) cell
%       - NSAMPLE : inteter
%       - 'SamplingMethodFunction' : function_handle
%       - 'SamplingMethodOptions' : (1,nOption) cell
%
%   Output:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - PADDINGSAMPLE : (nSampleExtra,nDesignVariable) double
%
%   See also sso_component_stochastic, candidate_space_sampling_all_feasible.
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

    parser = inputParser;
    parser.addParameter('SamplingMethodFunction',@sampling_latin_hypercube);
    parser.addParameter('SamplingMethodOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    nComponent = size(componentIndex,2);

    % create sampling box
    nDesignVariable = 0;
    for i=1:nComponent
        nDesignVariable = nDesignVariable + length(componentIndex{i});
    end
    samplingBox = nan(2,nDesignVariable);
    for i=1:nComponent
        samplingBox(:,componentIndex{i}) = candidateSpace(i).SamplingBox;
    end

    % start generating samples
    designSample = nan(nSample,size(samplingBox,2));
    paddingSampleComponent = cell(1,nComponent);
    nPadComponent = nan(1,nComponent);
    for i=1:nComponent
        designSampleComponent = [];
        paddingSampleComponent{i} = [];
        samplingBoxComponent = samplingBox(:,componentIndex{i});
        nGeneratedDesignSample = 0;

        % generate samples for each candidate space separately
        while (nGeneratedDesignSample <= nSample)
            currentSample = options.SamplingMethodFunction(samplingBoxComponent,nSample,options.SamplingMethodOptions{:});

            % see which designs are within current candidate spaces -> should be evaluated
            toBeEvaluated = candidateSpace(i).is_in_candidate_space(currentSample);
            nGeneratedDesignSample = nGeneratedDesignSample + sum(toBeEvaluated);

            % if we don't have enough, loop again;
            % if we do, get only those necessary and be done.
            if(nGeneratedDesignSample < nSample)
                if(sum(toBeEvaluated)>0)
                    designSampleComponent = [designSampleComponent;currentSample(toBeEvaluated,:)];
                    paddingSampleComponent{i} = [paddingSampleComponent{i};currentSample(~toBeEvaluated,:)];
                end

                continue;
            else
                % find how many extras must be removed and the main index
                nExtraSample = nGeneratedDesignSample - nSample;

                if(nExtraSample>0)
                    iExtra = find(toBeEvaluated,nExtraSample,'last');
                    iStopConcatenation = iExtra(1)-1;
                else
                    iStopConcatenation = length(toBeEvaluated);
                end

                % add everything before that
                designSampleComponent = [designSampleComponent;currentSample(toBeEvaluated(1:iStopConcatenation),:)];
                paddingSampleComponent{i} = [paddingSampleComponent{i};currentSample(~toBeEvaluated(1:iStopConcatenation),:)];
                break;
            end
        end
        
        designSample(:,componentIndex{i}) = designSampleComponent;
        nPadComponent(i) = size(paddingSampleComponent{i},1);
    end
    
    nPaddingBiggest = max(nPadComponent);
    nPadMissingOthers = nPaddingBiggest - nPadComponent;
    paddingSample = nan(nPaddingBiggest,size(samplingBox,2));
    for i=1:nComponent
        paddingAdditionalSample = [];
        needsAdditional = (nPadMissingOthers(i)>0);
        samplingBoxComponent = samplingBox(:,componentIndex{i});
        
        % sample inside bounding box with traditional methods
        paddingAdditionalSample = options.SamplingMethodFunction(samplingBoxComponent,nPadMissingOthers(i),options.SamplingMethodOptions{:});
        paddingSample(:,componentIndex{i}) = [paddingSampleComponent{i};paddingAdditionalSample];
    end
end