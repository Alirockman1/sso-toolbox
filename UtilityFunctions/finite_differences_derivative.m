function signalDerivative = finite_differences_derivative(signalValue,signalStep,derivativeOrder,nPoint)
%FINITE_DIFFERENCES_DERIVATIVE Derivative approximation using finite differences
%   FINITE_DIFFERENCES_DERIVATIVE numerically calculates the derivative of  
%   arbitrary order of an input signal using a finite differences approximation. 
%   The derivative is done by finite differences centered on each considered 
%   point; for left and right extremes, where centering is impossible, the same 
%   number of points is used, but the chosen samples are shifted to accomodate.
%
%   SIGNALDERIVATIVE = FINITE_DIFFERENCES_DERIVATIVE(SIGNALVALUE) calculates
%   the derivative SIGNALDERIVATIVE of the input signal SIGNALVALUE. It is
%   assumed the the interval between each sample (step size) is '1', that the 
%   first-order derivative is required, and that 7 points should be used to 
%   compute said derivative.
%
%   SIGNALDERIVATIVE = FINITE_DIFFERENCES_DERIVATIVE(SIGNALVALUE,SIGNALSTEP)
%   allows for the specification of the step-size used SIGNALSTEP. SIGNALSTEP
%   can also be given as an array of time with the same size as SIGNALVALUE,
%   in which case the mean difference between each sample will be used to
%   approximate the step size.
%
%   SIGNALDERIVATIVE = FINITE_DIFFERENCES_DERIVATIVE(SIGNALVALUE,SIGNALSTEP,
%   DERIVATIVEORDER) allows one to choose the derivative order DERIVATIVEORDER
%   of the output.
%
%   SIGNALDERIVATIVE = FINITE_DIFFERENCES_DERIVATIVE(SIGNALVALUE,SIGNALSTEP,
%   DERIVATIVEORDER,NPOINT) allows one to choose how many points will be used in
%   the finite differentes to approximate the derivative. Value is supposed to 
%   be an odd number, as derivatives are calculated centered around each point, 
%   which leads to symmetric relative indices, such as [-1,0,1] and 
%   [-2,-1,0,1,2], whose length must be odd-numbered.
%
%   Inputs:
%       - SIGNALVALUE : (1,nSample) double
%       - SIGNALSTEP : double OR (1,nSample) double
%       - DERIVATIVEORDER : integer
%       - NPOINT : integer
%
%   Output:
%       - SIGNALDERIVATIVE : (1,nSample) double
%
%   See also finite_differences_coefficients.
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

    if(nargin<4) 
        nPoint = [];
        if(nargin<3)
            derivativeOrder = [];
            if(nargin<2)
                signalStep = [];
            end
        end
    end
    signalStep = conditional_default_value_assignment(signalStep,[1,2]);
    derivativeOrder = conditional_default_value_assignment(derivativeOrder,1);
    nPoint = conditional_default_value_assignment(nPoint,7);

    % General information
    nSample = length(signalValue);
    stepSize = ternary(length(signalStep)>1,mean(diff(signalStep)),signalStep);
    denominator = stepSize^derivativeOrder;
    
    % Allocation
    signalDerivative = nan(size(signalValue));
    
    % Definition of left and right "extremes"
    leftLimit = floor(nPoint/2);
    rightLimit = nSample-floor(nPoint/2);
    
    % Derivatives on the left
    for i=1:leftLimit
        referenceLeft = 1-i;
        referenceRight = nPoint-i;
        sampledPoints = referenceLeft:1:referenceRight;
        coefficients = finite_differences_coefficients(sampledPoints,derivativeOrder);
        signalDerivative(i) = signalValue(i+sampledPoints)*coefficients/denominator;
    end
    
    % Derivatives on the right
    for i=nSample:-1:rightLimit
        referenceLeft = (nSample-i)-nPoint+1;
        referenceRight = (nSample-i);
        sampledPoints = referenceLeft:1:referenceRight;
        coefficients = finite_differences_coefficients(sampledPoints,derivativeOrder);
        signalDerivative(i) = signalValue(i+sampledPoints)*coefficients/denominator;
    end
    
    % Center derivatives
    nPointCenter = floor(nPoint/2);
    sampledPoints = -nPointCenter:1:nPointCenter;
    coefficients = finite_differences_coefficients(sampledPoints,derivativeOrder);
    for i=leftLimit+1:rightLimit-1
        signalDerivative(i) = signalValue(i+sampledPoints)*coefficients/denominator;
    end
    
    return;
end