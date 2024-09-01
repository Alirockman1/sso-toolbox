function coefficients = finite_differences_coefficients(sampledPoints,derivativeOrder)
%FINITE_DIFFERENCES_COEFFICIENTS Coefficients of each point for derivative
%   FINITE_DIFFERENCES_COEFFICIENTS finds the coefficients used for derivative
%   computation using finite differences for the given sample points.
%
%   COEFFICIENTS = FINITE_DIFFERENCES_COEFFICIENTS(SAMPLEDPOINTS,
%   DERIVATIVEORDER) uses the points with relative indexes SAMPLEDPOINTS to 
%   calculate the finite differences coefficients COEFFICIENTS for the 
%   derivative of order DERIVATIVEORDER. The relative indexes in SAMPLEDPOINTS 
%   refer to the point where the derivative is being calculated; for example, 
%   for an input of [-1,0,1] that means the points being used to compute a 
%   derivative are the values of the function at one point before, that point
%   itself, and one point after. COEFFICIENTS, then, would contain the 
%   respective coefficients to be used for those points (in that same example,
%   that could be [-1/2,0,1/2] for a first-order derivative, for example).
%
%   More information on how this is achieved can be found on:
%   https://web.media.mit.edu/~crtaylor/calculator.html
%
%   Inputs:
%       - SAMPLEDPOINTS : (1,nPoint) integer
%       - DERIVATIVEORDER : integer
%
%   Output:
%       - COEFFICIENTS : (nPoint,1) double
%
%   See also: finite_differences_derivative, factorial, mldivide.
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

    % sampled points must be integer; values must be unique and sorted
    if(any(mod(sampledPoints,1)~=0))
        coefficients = 0;
        return;
    end
    sampledPoints = unique(sampledPoints);
    nPoint = length(sampledPoints);
    
    % derivative order must be integer and positive
    if(mod(derivativeOrder,1)~=0 || derivativeOrder<1 || derivativeOrder>=nPoint)
        coefficients = 0;
        return;
    end
    
    % Build System Matrix Row-By-Row
    A = nan(nPoint,nPoint);
    for i=1:nPoint
        A(i,:) = sampledPoints.^(i-1);
    end
    
    % Build Right-Hand-Side Array
    b = zeros(nPoint,1);
    b(derivativeOrder+1) = factorial(derivativeOrder);
    
    % Solve
    coefficients = A\b;
end