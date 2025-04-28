function performanceMeasure = bead_slide_time_continuous(designSample,systemParameter)
%BEAD_SLIDE_TIME_CONTINUOUS Bottom-Up Mapping (Brachistochrone Problem)
%   BEAD_SLIDE_TIME_CONTINUOUS calculates the time a bead takes to slide down a
%   given curve. The curve is described by discrete points that are interpolated
%   using 'pchip' interpolation. The total time is calculated using trapezoidal
%   integration, assuming gravity of 9.80[m/s^2]. The optimal solution to this
%   problem is known as the brachistochrone curve.
%
%   PERFORMANCEMEASURE = BEAD_SLIDE_TIME_CONTINUOUS(DESIGNSAMPLE,
%   SYSTEMPARAMETER) receives the functionals used to determine height points of 
%   the curve in DESIGNSAMPLE and the system configuration in SYSTEMPARAMETER, 
%   returning the sliding time for each curve in PERFORMANCEMEASURE.
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,1) function_handle
%           -- function handle to the height points of the curve
%       - SYSTEMPARAMETER : (1,3) double
%           -- (1) : horizontal distance
%           -- (2) : vertical distance
%           -- (3) : number of interpolation points
%
%   Outputs:
%       - PERFORMANCEMEASURE : (nSample,1) double
%           -- sliding time for each curve
%
%   See also interp1, trapz, solve_analytical_brachistochrone.
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

    distanceX = systemParameter(1);
    distanceY = systemParameter(2);
    nInterpolationPoints = systemParameter(3);
    
    gravityConstant = 9.80;
    derivativeNeighborPoint = 2;
    
    nSample = size(designSample,1);
    performanceMeasure = nan(nSample,1);
    for i=1:nSample
        % create a base entry with the rough points given
        widthFine = linspace(0,distanceX,nInterpolationPoints);
        heightFine = designSample{i}(widthFine)-distanceY;

        % estimate the slope of the curve with finite differences
        slopeFine = finite_differences_derivative(heightFine,widthFine,1,derivativeNeighborPoint);
        
        % find the time derivative through the equations of motion
        timeDerivative = sqrt((1+slopeFine.^2)./(0-2*gravityConstant*heightFine));
        
        % First point is inf since heightFine(1)=0, approximate instead
        timeDerivative(1) = interp1(widthFine(2:end),timeDerivative(2:end),0,'pchip','extrap');
        
        % Integrate to find total time
        performanceMeasure(i) = trapz(widthFine,timeDerivative);
    end
end