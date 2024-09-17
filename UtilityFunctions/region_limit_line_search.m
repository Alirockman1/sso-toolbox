function [stepSizeLimitInside,stepSizeLimitOutside] = region_limit_line_search(regionCriterion,initialPoint,direction,designSpace,maxIter)
%REGION_LIMIT_LINE_SEARCH Find maximum step-sizes that remain inside region
%   REGION_LIMIT_LINE_SEARCH finds the maximum step-size that can be used for 
%   the given initial points and directions such that the new point is still 
%   inside the given region (considered individually the same as initial point)
%   and inside the design space. If no region criterion is given, only the 
%   design space is considered.
%   
%   STEPSIZELIMITINSIDE = REGION_LIMIT_LINE_SEARCH(REGIONCRITERION,INITIALPOINT,
%   DIRECTION,DESIGNSPACE) finds the maximum step sizes STEPSIZELIMITINSIDE 
%   such that, when starting from points INITIALPOINT and going in direction
%   DIRECTION, the new point is still inside the region defined by 
%   REGIONCRITERION (and also inside design space DESIGNSPACE). If 
%   REGIONCRITERION is empty, only the step size to remain inside the design
%   space is computed. It performs 20 iterations of both bracketing and 
%   bissectioning each to find said step sizes.
%
%   STEPSIZELIMITINSIDE = REGION_LIMIT_LINE_SEARCH(REGIONCRITERION,INITIALPOINT,
%   DIRECTION,DESIGNSPACE,MAXITER) allows one to specify the maximum number of
%   iterations MAXITER for the bracketing and bissectioning operations.
%
%   [STEPSIZELIMITINSIDE,STEPSIZELIMITOUTSIDE] = REGION_LIMIT_LINE_SEARCH(...)
%   also returns the minimum step size found to change from the inside of the 
%   region to the outside STEPSIZELIMITOUTSIDE. 
%
%   Input:
%       - REGIONCRITERION : function_handle
%       - INITIALPOINT : (nSample,nDesignVariable) double
%       - DIRECTION : (nSample,nDesignVariable) double
%       - DESIGNSPACE : (2,nDesignVariable) double
%       - MAXITER : integer
%
%   Output:
%       - STEPSIZELIMITINSIDE : (nSample,1) double
%       - STEPSIZELIMITOUTSIDE : (nSample,1) double
%
%   See also .
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

    if(nargin<5 || isempty(maxIter))
        maxIter = 20;
    end
    
    nSample = size(initialPoint,1);

    % initial step size -> see how much would be necessary to get to the edges 
    % of the design space assuming you are starting on the opposite one within 
    % maxIter iterations of bracketing
    stepSizeInitialDesignSpace = max(designSpace(2,:)-designSpace(1,:))./2^maxIter;
    stepSizeInitialDesignSpace = repmat(stepSizeInitialDesignSpace,nSample,1);

    % search for limit where design is still inside design space
    designSpaceCriterion = @(stepSize) [is_in_design_box(initialPoint + stepSize.*direction,designSpace)];
    [stepSizeLowerDesignSpace,stepSizeUpperDesignSpace] = bracketing_line_search(designSpaceCriterion,stepSizeInitialDesignSpace,maxIter);
    [stepSizeLowerDesignSpace,~] = bissection_line_search(designSpaceCriterion,stepSizeLowerDesignSpace,stepSizeUpperDesignSpace,maxIter);

    if(~isempty(regionCriterion))
        % given the maximum step size for design space, compute the initial one for
        % the given region in an analogous way
        stepSizeInitialRegion = stepSizeLowerDesignSpace./2^maxIter;
    
        % find bracketing limits for inside the given region
        regionStepSizeCriterium = @(stepSize) [regionCriterion(initialPoint + stepSize.*direction)];
        [stepSizeLowerRegion,stepSizeUpperRegion] = bracketing_line_search(regionStepSizeCriterium,stepSizeInitialRegion,maxIter); 
        [stepSizeLimitInside,stepSizeLimitOutside] = bissection_line_search(regionStepSizeCriterium,stepSizeLowerRegion,stepSizeUpperRegion,maxIter);
    else
        stepSizeLimitInside = stepSizeLowerDesignSpace;
        stepSizeLimitOutside = [];
    end
end

function [stepSizeLower,stepSizeUpper] = bracketing_line_search(regionCriterion,stepSize,maxIter)
	initialState = regionCriterion(0);
	previousStepSize = stepSize;
    previousState = regionCriterion(stepSize);
    
    nSample = size(stepSize,1);
    hasConverged = false(nSample,1);
    iterationBracketing = 1;
    while(~all(hasConverged) && (iterationBracketing<=maxIter))
    	% relocate previous result
        previousStepSize(~hasConverged) = stepSize(~hasConverged);

        % double/half step size
        mustIncrease = (initialState==previousState) & ~hasConverged;
        mustDecrease = (initialState~=previousState) & ~hasConverged;
        stepSize(mustIncrease) = 2*stepSize(mustIncrease);
        stepSize(mustDecrease) = 1/2*stepSize(mustDecrease);

        % evaluate
        currentState = regionCriterion(stepSize);

        % broad interval found if results between previous/current are different
        hasConverged = xor(previousState,currentState);
        
        iterationBracketing = iterationBracketing + 1;
    end

    stepSizeLower = min(previousStepSize,stepSize);
    stepSizeUpper = max(previousStepSize,stepSize);
end

function [stepSizeLower,stepSizeUpper] = bissection_line_search(regionCriterion,stepSizeLower,stepSizeUpper,maxIter)
	initialState = regionCriterion(0);

	% perform binary search to find approximate boundary
    for i=1:maxIter
    	% evaluate middle
        stepSize = stepSizeLower + (stepSizeUpper-stepSizeLower)./2;
        currentState = regionCriterion(stepSize);

        % move lower/upper boundary if positive/negative
        moveLower = (initialState==currentState);
        moveUpper = (initialState~=currentState);
        stepSizeLower(moveLower) = stepSize(moveLower);
        stepSizeUpper(moveUpper) = stepSize(moveUpper);
    end
end