function performanceMeasure = bead_slide_time(designSample,systemParameter)
%BEAD_SLIDE_TIME Bottom-up Mapping (Brachistochrone Problem)
%   BEAD_SLIDE_TIME returns the amount of time a bead takes to slide down a 
%   given curve. A discretized description of the curve is given and points in 
%   the middle are interpolated using a 'pchip' interpolation strategy. 
%   The total time is calculated using a simple integration strategy with 
%   trapezoidal approximations. Gravity has the assumed value of 9.80[m/s^2].
%   The optimum of this function is a brachistochrone.
%
%   PERFORMANCEMEASURE = BEAD_SLIDE_TIME(DESIGNSAMPLE,SYSTEMPARAMETER) 
%   receives the desired discretized height points of the curve at each division
%   in DESIGNSAMPLE, the initial horizontal/vertical distances and the number of
%   interpolation points in SYSTEMPARAMETER, and returns the sliding time for 
%   each curve / design sample point in PERFORMANCEMEASURE.
%
%   Example
%       Assuming a straight line between the starting and end point with
%       1m distance in both directions, 5 divisions and 10 
%       interpolation points:
%           distanceX = 1;
%           distanceY = 1;
%           nInterpolation = 10;
%           curve = linspace(distanceY,0,7); % straight line from start to finish
%           designSample = curve(2:6); % remove start/end point
%           systemParameter = [distanceX,distanceY,nInterpolation];
%           performanceMeasure = BEAD_SLIDE_TIME(designSample,systemParameter);
%       This will return the time it takes for a bead to roll this
%       straight line. 
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,nDivision) double
%       - SYSTEMPARAMETER : (1,3) double
%           -- (1) : horizontal distance
%           -- (2) : vertical distance
%           -- (3) : number of interpolation points
%
%   Outputs:
%       - PERFORMANCEMEASURE : (nSample,1) double
%
%   See also interp1, trapz, solve_analytical_brachistochrone.
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

    distanceX = systemParameter(1);
    distanceY = systemParameter(2);
    nInterpolationPoints = systemParameter(3);
    
    gravityConstant = 9.80;
    interpolationMethod = 'pchip';
    derivativeNeighborPoint = 2;
    
    nSample = size(designSample,1);
    performanceMeasure = nan(nSample,1);
    for i=1:nSample
        % create a base entry with the rough points given
        heightBase = [0,designSample(i,:)-distanceY,-distanceY]; % shift so start point is (0,0)
        widthBase = linspace(0,distanceX,size(heightBase,2));

        % create a fine mesh with the given interpolation points
        widthFine = linspace(0,distanceX,nInterpolationPoints);
        heightFine = interp1(widthBase,heightBase,widthFine,interpolationMethod);

        % estimate the slope of the curve with finite differences
        slopeFine = finite_differences_derivative(heightFine,widthFine,1,derivativeNeighborPoint);
        
        % find the time derivative through the equations of motion
        timeDerivative = sqrt((1+slopeFine.^2)./(0-2*gravityConstant*heightFine));
        
        % First point is inf since heightFine(1)=0, approximate instead
        timeDerivative(1) = timeDerivative(2);
        
        % Integrate to find total time
        performanceMeasure(i) = trapz(widthFine,timeDerivative);
    end
end