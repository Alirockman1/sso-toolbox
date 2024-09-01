function [designSample,paddingSample] = candidate_space_sampling_all_feasible(candidateSpace,component,nSample,varargin)
%CANDIDATE_SPACE_SAMPLING_ALL_FEASIBLE Sampling inside candidate spaces
%   CANDIDATE_SPACE_SAMPLING_ALL_FEASIBLE produces design sample points that
%   are inside all of the candidates spaces; it does so by generating a sample
%   in the complete space and then checking which design points are 
%   simultaneously inside all of the candidate spaces. While not enough samples
%   have been generated, the process repeats.
%
%   DESIGNSAMPLE = CANDIDATE_SPACE_SAMPLING_ALL_FEASIBLE(CANDIDATESPACE,
%   COMPONENT,NSAMPLE) generates NSAMPLE design sample points that are inside 
%   the candidate spaces CANDIDATESPACE, with components COMPONENT, returning
%   said sample in DESIGNSAMPLE.
%
%   DESIGNSAMPLE = CANDIDATE_SPACE_SAMPLING_ALL_FEASIBLE(...NAME,VALUE,...) 
%   allows for the specification of additional options. These are:
%       - 'SamplingMethodFunction' : base sampling method to be used. Default:
%       @sampling_latin_hypercube.
%       - 'SamplingMethodOptions' : extra options for the base sampling method.
%       Default is empty.
%
%   [DESIGNSAMPLE,PADDINGSAMPLE] = CANDIDATE_SPACE_SAMPLING_ALL_FEASIBLE(...)
%   additionally returns extra samples generated PADDINGSAMPLE, which are not
%   inside at least one of the candidate spaces.
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
%   See also sso_component_stochastic, 
%   candidate_space_sampling_individual_feasible.
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

    nComponent = size(component,2);

    % create sampling box
    nDesignVariable = 0;
    for i=1:nComponent
        nDesignVariable = nDesignVariable + length(component{i});
    end
    samplingBox = nan(2,nDesignVariable);
    for i=1:nComponent
        samplingBox(:,components{i}) = candidateSpace(i).SamplingBox;
    end

    % start generating samples
    designSample = [];
    paddingSample = [];
    nGeneratedDesignSample = 0;
    while (nGeneratedDesignSample <= nSample)
        % sample inside bounding box with traditional methods
        initialsample = options.SamplingMethodFunction(samplingBox,nSample,options.SamplingMethodOptions{:});
        
        % see which designs are within all candidate spaces -> should be evaluated
        isInsideCandidateSpace = false(nSample,nComponent);
        for i=1:nComponent
            isInsideCandidateSpace(:,i) = candidateSpace(i).is_in_candidate_space(initialsample(:,component{i}));
        end
        toBeEvaluated = all(isInsideCandidateSpace,2);
        nGeneratedDesignSample = nGeneratedDesignSample + sum(toBeEvaluated);
        
        % if we don't have enough, loop again;
        % if we do, get only those necessary and be done.
        if(nGeneratedDesignSample < nSample)
            if(sum(toBeEvaluated)>0)
                designSample = [designSample;initialsample(toBeEvaluated,:)];
                paddingSample = [paddingSample;initialsample(~toBeEvaluated,:)];
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
            designSample = [designSample;initialsample(toBeEvaluated(1:iStopConcatenation),:)];
            paddingSample = [paddingSample;initialsample(~toBeEvaluated(1:iStopConcatenation),:)];
            break;
        end
    end
end