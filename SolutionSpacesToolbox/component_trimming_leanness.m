function activeComponent = component_trimming_leanness(designSample,labelKeep,trimmingOrder,component,varargin)
%COMPONENT_TRIMMING_LEANNESS Apply the leanness condition to the current problem
%   COMPONENT_TRIMMING_LEANNESS trims out all design points that can be removed 
%   without trimming out designs that should be kept. This is mainly used to
%   remove design points labeled as useless without removing any designs
%   labeled as useful, as required by the leanness condition of requirement
%   spaces. 
%
%   ACTIVECOMPONENT = COMPONENT_TRIMMING_LEANNESS(DESIGNSAMPLE,LABELKEEP,
%   TRIMMINGORDER,COMPONENT) uses the design sample points in DESIGNSAMPLE 
%   and trims out (if possible) the designs specified in TRIMMINGORDER while 
%   keeping the designs with that shouldn't be removed (as specified in 
%   LABELKEEP). The result is returned in ACTIVECOMPONENT for each component,
%   with labels as to whether each design point is inside (true) or outside 
%   (false) the component space.
%
%   ACTIVECOMPONENT = COMPONENT_TRIMMING_LEANNESS(DESIGNSAMPLE,LABELKEEP,
%   TRIMMINGORDER,COMPONENT,ACTIVECOMPONENT) allows one to specify an initial
%   state for ACTIVECOMPONENT, so not all designs are assumed to be inside
%   at the start.
%
%   ACTIVECOMPONENT = COMPONENT_TRIMMING_LEANNESS(...NAME,VALUE,...) allows the 
%   specification of name-value pair arguments. These can be:
%       - 'TrimmingMethodFunction' : function handle of the method used for the 
%       trimming operation. Default is 
%       'component_trimming_method_planar_trimming'.
%       - 'TrimmingMethodOptions' : any extra options for the trimming method.
%       Default is empty.
%
%   Input:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - LABELKEEP : (nSample,1) logical
%       - TRIMMINGORDER : (nExclude,1) integer
%       - COMPONENT : (1,nComponent) cell
%       - ACTIVECOMPONENT : (nSample,nComponent) logical
%       - 'TrimmingMethodFunction' : function_handle
%       - 'TrimmingMethodOptions' : (1,nOption) cell
%
%   Output:
%       - ACTIVECOMPONENT : (nSample,nComponent) logical
%
%   See also component_trimming_operation.
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
    parser.KeepUnmatched = true;
    parser.addRequired('designSample',@(x)isnumeric(x));
    parser.addRequired('labelKeep',@(x)islogical(x));
    parser.addRequired('trimmingOrder',@(x)isnumeric(x));
    parser.addRequired('component',@(x)iscell(x)&&(size(x,1)==1));
    parser.addOptional('activeComponent',[],@(x)islogical(x)||isempty(x));
    parser.addParameter('TrimmingMethodFunction',@component_trimming_method_planar_trimming,@(x)isa(x,'function_handle'));
    parser.addParameter('TrimmingMethodOptions',{},@(x)(iscell(x)));
    parser.parse(designSample,labelKeep,trimmingOrder,component,varargin{:});

    % unwrap
    activeComponent = parser.Results.activeComponent;
    trimmingMethodFunction = parser.Results.TrimmingMethodFunction;
    trimmingMethodOptions = parser.Results.TrimmingMethodOptions;

    nSample = size(designSample,1);
    nExclude = size(trimmingOrder,1);
    nComponent = size(component,2);
    activeComponent = conditional_default_value_assignment(activeComponent,true(nSample,nComponent));
    activeAll = all(activeComponent,2);
    activeKeep = activeAll & labelKeep;
    
    for i=1:nExclude
        iExclude = trimmingOrder(i);
        if(~activeAll(iExclude))
            continue;
        end

        for j=1:nComponent
            activeComponentDesign = activeComponent(:,j);

            designSampleComponent = designSample(activeComponentDesign,component{j});
            iRemovalComponent = convert_index_base(activeComponentDesign,iExclude,'forward');
            activeKeepComponent = activeKeep(activeComponentDesign);

            removalCandidateComponent = trimmingMethodFunction(...
                designSampleComponent,iRemovalComponent,activeKeepComponent,trimmingMethodOptions{:});

            removalCandidate = false(nSample,size(removalCandidateComponent,2));
            removalCandidate(activeComponentDesign,:) = removalCandidateComponent;

            % see if any elimination is possible without eliminating designs that must be kept
            unwantedRemoval = removalCandidate & activeKeep;
            canBeRemoved = ~any(unwantedRemoval,1);
            trimRemoval = any(removalCandidate(:,canBeRemoved),2);

            activeComponent(trimRemoval,j) = false;
            activeAll(trimRemoval) = false;
        end
    end
end