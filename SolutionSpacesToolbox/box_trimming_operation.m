function [candidateBoxTrimmed,measureTrimmed] = box_trimming_operation(designSample,isViable,trimmingOrder,varargin)
%BOX_TRIMMING_OPERATION performs the box trimming operation for solution spaces
%   BOX_TRIMMING_OPERATION performs the trimming operation of the design box 
%   candidate, in such a way that the final result only includes good 
%   (acceptable) designs and has the biggest possible measure.
%   
%   CANDIDATEBOXTRIMMED = BOX_TRIMMING_OPERATION(DESIGNSAMPLE,ISVIABLE,
%   TRIMMINGORDER) receives the design sample points in DESIGNSAMPLE, a label 
%   indicating if they should be kept ISVIABLE, and the order of designs to 
%   remove TRIMMINGORDER, and returns the optimal design box with only 
%   acceptable designs in CANDIDATEBOXTRIMMED.
%
%   CANDIDATEBOXTRIMMED = BOX_TRIMMING_OPERATION(DESIGNSAMPLE,ISVIABLE,
%   TRIMMINGORDER,CANDIDATEBOX) allows the specification of an initial candidate
%   box. If not given, a bounding box is used around DESIGNSAMPLE.
%
%   CANDIDATEBOXTRIMMED = BOX_TRIMMING_OPERATION(...,NAME,VALUE,...) allows for 
%   the choice of specific optimization parameters, as opposed to the use of
%   the default ones. This options available are:
%       - 'MeasureFunction' : function to calculate the measure of each 
%       candidate box. Default: @box_measure_volume.
%       - 'MeasureOptions' : options to be used in the measure function. Default
%       is empty.
%       - 'Slack' : a value between 0 and 1 which specifies where the boundary
%       should be relocated to when performing the trimming operation: for 0,
%       the boundary of the design box gets relocated to the closest acceptable
%       design; for 1, the boundary gets relocated to the unnaceptable design
%       being removed; and proportionally in-between for values between 0 and 1.
%       - 'PassesCriterion' : criterion for what passes to consider in terms
%       of the good anchors considered; in each pass, the chosen viable point
%       cannot be removed in the operation, and at the end of the pass, the 
%       result from that gets compared to the previous optimal trimming. The 
%       following options are available:
%           -- 'single' : no viable anchors are considered, and only a single 
%           pass is done
%           -- 'reduced' : all viable points must be included inside  at least
%           one solution.
%           -- 'full' : one individual solution is generated for each viable 
%           point.
%
%   [CANDIDATEBOXTRIMMED,MEASURETRIMMED] = BOX_TRIMMING_OPERATION(...) also 
%   returns the calcuated measure of the optimal trimmed box.
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - ISVIABLE : (nSample,1) logical
%       - TRIMMMINGORDER : 
%       - CANDIDATEBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - 'MeasureFunction' : function_handle with interface 
%       'measure = f(candidateBox,weight,nSample,nUse)'
%       - 'MeasureOptions' : (1,nOptions) cell
%       - 'Slack' : double
%       - 'PassesCriterion' : char OR string
%
%   Outputs:
%       - CANDIDATEBOXTRIMMED : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - MEASURETRIMMED : double
%
%   See also sso_box_stochastic, box_measure_volume.
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

    % parse inputs
    parser = inputParser;
    parser.addRequired('designSample',@(x)isnumeric(x));
    parser.addRequired('isViable',@(x)islogical(x));
    parser.addRequired('trimmingOrder',@(x)isnumeric(x));
    parser.addOptional('CandidateBox',[],@(x)isnumeric(x)&&(size(x,1)==2 || isempty(x)));
    parser.addParameter('MeasureFunction',@box_measure_volume,@(x)isa(x,'function_handle'));
    parser.addParameter('MeasureOptions',{},@(x)iscell(x));
    parser.addParameter('Slack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
    parser.addParameter('PassesCriterion','full',@(x)any(strcmpi(x,{'single','reduced','full'})));
    
    parser.parse(designSample,isViable,trimmingOrder,varargin{:});
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

    isInsideBoxInitial = is_in_design_box(designSample,candidateBoxInitial);
    iAnchorViable = [];
    canBeAnchorViable = isViable & isInsideBoxInitial;
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

                        closestViable = [];
                        if(options.Slack<1)
                            isInsideBox = is_in_design_box(designSample,candidateBox);
                            if(l==1) % lower boundary
                                remainRegion = (designSample(:,k)>=designSample(iExclude,k));
                                closestViable = min(designSample(isInsideBox&remainRegion&isViable,k));
                            else % upper boundary
                                remainRegion = (designSample(:,k)<=designSample(iExclude,k));
                                closestViable = max(designSample(isInsideBox&remainRegion&isViable,k));
                            end
                        end

                        if(isempty(closestViable))
                            closestViable = designSample(iExclude,k);
                        end

                        % move lower/upper bound
                        candidateBox(l,k) = designSample(iExclude,k)*options.Slack + closestViable*(1-options.Slack);
                        isInsideBox = is_in_design_box(designSample,candidateBox);

                        % if anchor is not included, do not use; otherwise, compute measure
                        if(~isInsideBox(iAnchorViable))
                            measureBox = -inf;
                        else
                            isInsideBoxViable = isInsideBox & isViable;
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

            % check final measure
            isInsideBox = is_in_design_box(designSample,candidateBoxCurrent);
            isInsideBoxViable = isInsideBox & isViable;
            fractionViable = sum(isInsideBoxViable)/sum(isInsideBox);
            measureBox = options.MeasureFunction(candidateBoxCurrent,fractionViable,options.MeasureOptions{:});

            % check final measure
            if(measureBox>measureTrimmed)
                candidateBoxTrimmed = candidateBoxCurrent;
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