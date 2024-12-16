function trimmedCandidateSpace = component_trimming_operation(designSample,labelViable,trimmingOrder,componentIndex,candidateSpace,varargin)
%COMPONENT_TRIMMING_OPERATION Component SSO Removal of Unnaceptable Designs
%   COMPONENT_TRIMMING_OPERATION is the main trimming operation for component 
%   solution space optimization. In this, all design sample points marked for 
%   removal are excluded from the component space, with the trimming method
%   specified. This procedure can accept multiple trimming orders, and tries 
%   them procedurally. Additionally, if set accordingly, it can also produce 
%   results which always keep a particular design one at a time, and then choose
%   the best overall result from those.
%
%   ACTIVECOMPONENT = COMPONENT_TRIMMING_OPERATION(DESIGNSAMPLE,LABELVIABLE,
%   TRIMMINGORDER,COMPONENT) uses the design sample points in DESIGNSAMPLE 
%   and trims out the designs specified in TRIMMINGORDER, taking into 
%   consideration the designs with that shouldn't be removed (as specified in 
%   LABELVIABLE). The result is returned in ACTIVECOMPONENT for each component,
%   with labels as to whether each design point is inside (true) or outside 
%   (false) the component space.
%
%   ACTIVECOMPONENT = COMPONENT_TRIMMING_OPERATION(DESIGNSAMPLE,LABELVIABLE,
%   TRIMMINGORDER,COMPONENT,ACTIVECOMPONENT) allows one to specify an initial
%   state for ACTIVECOMPONENT, so not all designs are assumed to be inside
%   at the start.
%
%   ACTIVECOMPONENT = COMPONENT_TRIMMING_OPERATION(...NAME,VALUE,...) allows the 
%   specification of name-value pair arguments. These can be:
%       - 'TrimmingMethodFunction' : function handle of the method used for the 
%       trimming operation. Default is 
%       'component_trimming_method_planar_trimming'.
%       - 'TrimmingMethodOptions' : any extra options for the trimming method.
%       Default is empty.
%       - 'TrimmingCostFunction' : function handle of the method used to 
%       determine the cost of each removal candidate for the considered 
%       component. Default is 'component_trimming_cost'.
%       - 'TrimmingCostOptions' : any extra options for the the calculation of
%       cost of each removal candidate. Default is empty.
%       - 'TrimmingComponentChoiceFunction' : function handle of the method used  
%       to determine the component trim to be used given the component costs. 
%       Default is 'component_trimming_choice'.
%       - 'TrimmingComponentChoiceOptions' : any extra options for the choice
%       of component trim to perform. Default is empty.
%       - 'PassesCriterion' : choice regarding whether or not the operation 
%       should try using anchor points (points that must be kept and cannot
%       be trimmed). Possible values are:
%           -- 'single' : no anchor points are considered; only a single pass
%           of the trimming operation is done.
%           -- 'reduced' : after each pass is done, all designs included
%           in the final result are marked as having been used as anchors. 
%           The process repeats itself until all design points have been 
%           marked. The best overall is then selected.
%           -- 'full' : process repeats once for each design that should be 
%           kept, considering that as anchor. The best overall is then selected.
%
%   Input:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - LABELVIABLE : (nSample,1) logical
%       - TRIMMINGORDER : (nExlude,nOrder) integer
%       - COMPONENT : (1,nComponent) cell
%       - ACTIVECOMPONENT : (nSample,nComponent) logical
%       - 'TrimmingMethodFunction' : function_handle
%       - 'TrimmingMethodOptions' : (1,nOptionsMethod) cell
%       - 'TrimmingCostFunction' : function_handle
%       - 'TrimmingCostOptions' : (1,nOptionsCost) cell
%       - 'TrimmingComponentChoiceFunction' : function_handle
%       - 'TrimmingComponentChoiceOptions' : (1,nOptionsChoice) cell
%       - 'PassesCriterion' : char OR string
%
%   Output:
%       - ACTIVECOMPONENT : (nSample,nComponent) logical
%
%   See also component_trimming_leanness, sso_component_stochastic, 
%   trimming_order.
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
    parser.addRequired('componentIndex',@(x)iscell(x)&&(size(x,1)==1));
    parser.addRequired('candidateSpace',@(x)isa(x,'CandidateSpaceBase'));
    parser.addParameter('TrimmingMethodFunction',@component_trimming_method_planar_trimming,@(x)isa(x,'function_handle'));
    parser.addParameter('TrimmingMethodOptions',{},@(x)(iscell(x)));
    parser.addParameter('TrimmingCostFunction',@component_trimming_cost,@(x)isa(x,'function_handle'));
    parser.addParameter('TrimmingCostOptions',{},@(x)(iscell(x)));
    parser.addParameter('TrimmingComponentChoiceFunction',@component_trimming_choice,@(x)isa(x,'function_handle'));
    parser.addParameter('TrimmingComponentChoiceOptions',{},@(x)(iscell(x)));
    parser.addParameter('PassesCriterion','reduced',@(x)any(strcmpi(x,{'single','reduced','full'})));
    parser.parse(designSample,labelViable,trimmingOrder,componentIndex,candidateSpace,varargin{:});

    % unwrap
    trimmingMethodFunction = parser.Results.TrimmingMethodFunction;
    trimmingMethodOptions = parser.Results.TrimmingMethodOptions;
    trimmingCostFunction = parser.Results.TrimmingCostFunction;
    trimmingCostOptions = parser.Results.TrimmingCostOptions;
    trimmingComponentChoiceFunction = parser.Results.TrimmingComponentChoiceFunction;
    trimmingComponentChoiceOptions = parser.Results.TrimmingComponentChoiceOptions;
    passesCriterion = parser.Results.PassesCriterion;

    nSample = size(designSample,1);
    nExclude = size(trimmingOrder,1);
    nTrimmingOrder = size(trimmingOrder,2);
    nComponent = size(componentIndex,2);

    % find which samples are already inside/outside each candidate space
    activeComponentInitial = true(size(designSample,1),nComponent);
    for i=1:size(componentIndex,2)
        activeComponentInitial(:,i) = candidateSpace(i).is_in_candidate_space(designSample(:,componentIndex{i}));
    end
    activeAllInitial = all(activeComponentInitial,2);
    activeKeepInitial = activeAllInitial & labelViable;

    totalCostMinimum = inf;
    optimalActiveComponent = activeComponentInitial;
    optimalTrimmingInformation = cell(1,nComponent);
    
    iAnchorViable = [];
    canBeAnchorViable = labelViable & activeAllInitial;
    hasBeenAnchorViable = false(1,sum(canBeAnchorViable));
    
    iCurrentPass = 0;
    doneWithPasses = isempty(trimmingOrder);
    while(~doneWithPasses)
        for i=1:nTrimmingOrder
            activeComponentPass = activeComponentInitial;
            activeAllPass = activeAllInitial;
            trimTotalPass = false(nSample,1);
            trimmingInformation = cell(1,nComponent);

            for j=1:nExclude
                iExclude = trimmingOrder(j,i);
                if(~activeAllPass(iExclude))
                    continue;
                end
                
                activeKeep = activeAllPass & labelViable;
                componentCost = nan(1,nComponent);
                componentRemoval = false(nSample,nComponent);
                for k=1:nComponent
                    activeComponentDesign = activeComponentPass(:,k);

                    designSampleComponent = designSample(activeComponentDesign,componentIndex{k});
                    iRemovalComponent = convert_index_base(activeComponentDesign,iExclude,'forward');
                    activeKeepComponent = activeKeep(activeComponentDesign);

                    [removalCandidateComponent,trimmingInformationCandidateComponent] = trimmingMethodFunction(...
                        designSampleComponent,iRemovalComponent,activeKeepComponent,trimmingMethodOptions{:});

                    removalCandidate = false(nSample,size(removalCandidateComponent,2));
                    removalCandidate(activeComponentDesign,:) = removalCandidateComponent;

                    violateAnchor = (removalCandidate(iAnchorViable,:)==true);
                    [removalCost,iMinimumCostRemoval] = trimmingCostFunction(...
                        designSample,activeKeep,removalCandidate,violateAnchor,trimmingCostOptions{:});
                    
                    componentCost(k) = removalCost;
                    componentRemoval(:,k) = removalCandidate(:,iMinimumCostRemoval);
                    trimmingInformationComponent(k) = trimmingInformationCandidateComponent(iMinimumCostRemoval);
                end

                % decide between components
                iComponentTrim = trimmingComponentChoiceFunction(...
                    componentCost,componentIndex,activeKeep,componentRemoval,trimmingComponentChoiceOptions{:});
                
                % Eliminate those designs
                trimRemoval = componentRemoval(:,iComponentTrim);
                activeComponentPass(trimRemoval,iComponentTrim) = false;
                activeAllPass(trimRemoval) = false;
                trimTotalPass(trimRemoval) = true;
                trimmingInformation{iComponentTrim} = [trimmingInformation{iComponentTrim};trimmingInformationComponent(iComponentTrim)];
            end

            % check total cost
            % FIX perhaps use measure here instead of total cost
            totalCostPass = trimmingCostFunction(designSample,activeKeepInitial,trimTotalPass,[],trimmingCostOptions{:});
            if(totalCostPass<totalCostMinimum)
                optimalActiveComponent = activeComponentPass;
                optimalTrimmingInformation = trimmingInformation;
                totalCostMinimum = totalCostPass;
            end
        end

        % check for convergence / update for next pass
        activeComponentViable = all(activeComponentPass,2) & labelViable;
        if(strcmpi(passesCriterion,'single'))
            doneWithPasses = true;
        else
            if(strcmpi(passesCriterion,'reduced'))
                includedInFinalResult = convert_index_base(canBeAnchorViable,activeComponentViable,'forward');
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

    for i=1:size(componentIndex,2)
        % Eliminate from sampling designs that were taken out by other components
        designSampleComponent = designSample(:,componentIndex{i});
        isInsideComponent = optimalActiveComponent(:,i);
        trimmingInformation = optimalTrimmingInformation{i};
        trimmedCandidateSpace(i) = candidateSpace(i).update_candidate_space(designSampleComponent,isInsideComponent,trimmingInformation);
    end
end