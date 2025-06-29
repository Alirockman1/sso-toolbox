function [candidateBoxTrimmed,measureTrimmed] = box_trimming_operation(designSample,labelViable,trimmingOrder,varargin)
%box_trimming_classic performs the box trimming operation for solution spaces
%   box_trimming_classic performs the trimming operation of the design box 
%   candidate, in such a way that the final result only includes good 
%   (acceptable) designs and has the biggest possible measure.
%   
%   CANDIDATEBOXTRIMMED = box_trimming_classic(CANDIDATEBOX,DESIGNSAMPLE,
%   LABELACCEPT,LABELUSE,SCORE) receives the current candidate design box in
%   CANDIDATEBOX, the design samples in DESIGNSAMPLE, their respective 
%   requirement spaces labels of acceptance LABELACCEPT and usefulness LABELUSE,
%   and their design score in SCORE, and returns the optimal design box with
%   only acceptable designs in CANDIDATEBOXTRIMMED. For LABELUSE and 
%   LABELACCEPT, 'true' means the design is useful/acceptable, and 'false' means
%   the opposite. SCORE is used to define the order of trimming.
%
%   CANDIDATEBOXTRIMMED = box_trimming_classic(CANDIDATEBOX,DESIGNSAMPLE,
%   LABELACCEPT,LABELUSE,SCORE,OPTIONS) allows for the choice of specific
%   optimization parameters as described in OPTIONS, as opposed to the use of
%   the default ones. This options available are:
%       - 'BoxMeasureFunction': function to calculate the measure of each 
%       candidate box.
%       - 'Weight': weight for each dimension, where dimensions with higher 
%       weights are given bigger priority.
%       - 'ApplyLeannss': a boolean flag to determine if the leanness condition
%       should be applied at each trimming step or not.
%       - 'MustIncludeReferenceDesign': a boolean flag to determine if the
%       reference design should always be a part of the candidate design box.
%       If 'true', design boxes that do not contain the reference design will be
%       given '-inf' measure, and therefore not be considered.
%       - 'ReferenceDesign': reference design, which may be required to be 
%       inside the solution box.
%       - 'Slack': a value between 0 and 1 which specifies where the boundary
%       should be relocated to when performing the trimming operation: for 0,
%       the boundary of the design box gets relocated to the closest acceptable
%       design; for 1, the boundary gets relocated to the unnaceptable design
%       being removed; and proportionally in-between for values between 0 and 1.
%       - 'Order': the criterium for the order of the trimming operation, 
%       specifically the order which unnacceptable designs must be removed.
%       Available options are 'score-low-to-high', 'score-high-to-low' and
%       'random': in the first, designs closest to the boundary are removed 
%       first; in the second, the designs furthest from the boundary are removed
%       first; and in the last, the order is randomized.
%
%   [CANDIDATEBOXTRIMMED,MEASURETRIMMED] = box_trimming_classic(...) also 
%   returns the calcuated measure of the optimal trimmed box.
%
%   Inputs:
%       - CANDIDATEBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - LABELACCEPT : (nSample,1) logical
%       - LABELUSE : (nSample,1) logical 
%       - SCORE : (nSample,1) double
%       - OPTIONS : structure OR name-value pair arguments, optional
%           -- BoxMeasureFunction : function handle with interface 
%           'measure = f(candidateBox,weight,nSample,nUse)'
%           -- weight : (1,nDesignVariable) double
%           -- ApplyLeanness : logical
%           -- MustIncludeReferenceDesign : logical
%           -- ReferenceDesign : (1,nDesignVariable) double
%           -- Slack : double
%           -- Order : string
%
%   Outputs:
%       - CANDIDATEBOXTRIMMED : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - MEASURETRIMMED : double
%
%   See also sso_box_stochastic, box_measure_volume.
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

    % parse inputs
    parser = inputParser;
    parser.addRequired('designSample',@(x)isnumeric(x));
    parser.addRequired('labelViable',@(x)islogical(x));
    parser.addRequired('trimmingOrder',@(x)isnumeric(x));
    parser.addOptional('CandidateBox',[],@(x)isnumeric(x)&&(size(x,1)==2 || isempty(x)));
    parser.addParameter('MeasureFunction',@box_measure_volume,@(x)isa(x,'function_handle'));
    parser.addParameter('MeasureOptions',{},@(x)iscell(x));
    parser.addParameter('Slack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
    parser.addParameter('PassesCriterion','full',@(x)any(strcmpi(x,{'single','reduced','full'})));
    
    parser.parse(designSample,labelViable,trimmingOrder,varargin{:});
    options = parser.Results;

    % assign conditional default values
    candidateBoxInitial = conditional_default_value_assignment(options.CandidateBox,design_bounding_box(designSample));
    
    % initialization
    nSample = size(designSample,1);
    nDesignVariable = size(designSample,2);
    nExclude = size(trimmingOrder,1);
    nTrimmingOrder = size(trimmingOrder,2);

    
    candidateBoxTrimmed = candidateBoxInitial;
    measureTrimmed = -inf;
    % isInsideBox = is_in_design_box(designSample,candidateBoxTrimmed);
    % isInsideBoxViable = labelViable & isInsideBox;
    % measureTrimmed = options.MeasureFunction(candidateBoxTrimmed,sum(isInsideBoxViable),options.MeasureOptions{:});

    isInsideBoxInitial = is_in_design_box(designSample,candidateBoxInitial);
    iAnchorViable = [];
    canBeAnchorViable = labelViable & isInsideBoxInitial;
    hasBeenAnchorViable = false(1,sum(canBeAnchorViable));

    iCurrentPass = 0;
    doneWithPasses = isempty(trimmingOrder);
    while(~doneWithPasses)
        for i=1:nTrimmingOrder
            candidateBoxCurrent = candidateBoxInitial;
            isInsideBoxCurrent = isInsideBoxInitial;

            for j=1:nExclude
                iExclude = trimmingOrder(j,i);
                if(~isInsideBoxCurrent(iExclude))
                    continue;
                end

                candidateBoxBest = [];
                isInsideBoxBest = [];
                measureBoxBest = -inf;
                
                for k=1:nDesignVariable
                    for l=[1,2]
                        candidateBox = candidateBoxCurrent;

                        % move lower/upper bound
                        candidateBox(l,k) = designSample(iExclude,k);
                        isInsideBox = is_in_design_box(designSample,candidateBox);

                        % if anchor is not included, do not use; otherwise, compute measure
                        if(~isInsideBox(iAnchorViable))
                            measureBox = -inf;
                        else
                            isInsideBoxViable = isInsideBox & labelViable;
                            fractionViable = sum(isInsideBoxViable)/sum(isInsideBox);
                            measureBox = options.MeasureFunction(candidateBox,fractionViable,options.MeasureOptions{:});
                        end

                        if(measureBox > measureBoxBest)
                            candidateBoxBest = candidateBox;
                            isInsideBoxBest = isInsideBox;
                            measureBoxBest = measureBox;
                        end
                    end
                end

                % update 
                candidateBoxCurrent = candidateBoxBest;
                isInsideBoxCurrent = isInsideBoxBest;
            end

            % use slack
            strictBox = design_bounding_box(designSample,isInsideBoxCurrent & labelViable);
            candidateBox = candidateBoxCurrent*options.Slack + strictBox*(1-options.Slack);
            isInsideBox = is_in_design_box(designSample,candidateBox);
            isInsideBoxViable = isInsideBox & labelViable;
            fractionViable = sum(isInsideBoxViable)/sum(isInsideBox);
            measureBox = options.MeasureFunction(candidateBox,fractionViable,options.MeasureOptions{:});

            % check final measure
            if(measureBox>measureTrimmed)
                candidateBoxTrimmed = candidateBox;
                measureTrimmed = measureBox;
            end
        end

        % check for convergence / update for next pass
        if(strcmpi(options.PassesCriterion,'single'))
            doneWithPasses = true;
        else
            if(strcmpi(options.PassesCriterion,'reduced'))
                includedInFinalResult = convert_index_base(canBeAnchorViable,isInsideBoxViable,'forward');
                hasBeenAnchorViable(includedInFinalResult) = true;

                currentAnchor = convert_index_base(canBeAnchorViable,iAnchorViable,'forward');
                hasBeenAnchorViable(currentAnchor) = true; % even if it was removed (no choice)
            else %if(strcmpi(options.PassesCriterion,'full'))
                if(iCurrentPass~=0)
                    hasBeenAnchorViable(iCurrentPass) = true;
                end
            end

            notYetAnchor = find(~hasBeenAnchorViable,1,'first');
            if(isempty(notYetAnchor))
                doneWithPasses = true;
            else
                iAnchorViable = convert_index_base(canBeAnchorViable,notYetAnchor,'backward');
            end
        end
        iCurrentPass = iCurrentPass + 1;
    end
end