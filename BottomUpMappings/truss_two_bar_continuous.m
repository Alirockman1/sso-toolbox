function performanceMeasure = truss_two_bar_continuous(designSample,systemParameter)
%TRUSS_TWO_BAR_CONTINUOUS Bottom-up Mapping (Truss Stiffness Problem)
%   TRUSS_TWO_BAR_CONTINUOUS calculates the performance measure for a truss
%   with two bars connected in parallel. The function computes the total
%   stiffness of the system by interpolating variable material properties
%   (Young's modulus and cross-sectional area) along the length of each bar.
%   The stiffness is calculated using numerical integration with trapezoidal
%   approximations.
%
%   PERFORMANCEMEASURE = TRUSS_TWO_BAR_CONTINUOUS(DESIGNSAMPLE,SYSTEMPARAMETER)
%   receives the material properties for both bars at discrete points in
%   DESIGNSAMPLE, the lengths of both bars and the number of interpolation
%   points in SYSTEMPARAMETER, and returns the total stiffness for each
%   design sample in PERFORMANCEMEASURE.
%
%   Example
%       Assuming two bars with lengths 1m and 1.5m, using 5 divisions and
%       10 interpolation points:
%           length1 = 1;
%           length2 = 1.5;
%           nInterpolation = 10;
%           nDivisions = 5;
%           % Create sample with constant properties for both bars
%           youngsModulus1 = ones(1,nDivisions)*200e9; % 200 GPa
%           area1 = ones(1,nDivisions)*0.001; % 0.001 m²
%           youngsModulus2 = ones(1,nDivisions)*200e9; % 200 GPa
%           area2 = ones(1,nDivisions)*0.001; % 0.001 m²
%           designSample = [youngsModulus1, area1, youngsModulus2, area2];
%           systemParameter = [length1, length2, nInterpolation];
%           performanceMeasure = truss_two_bar_continuous(designSample,systemParameter);
%       This will return the total stiffness of the two-bar truss system.
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,4*nDivisions) double
%           -- (1:nDivisions) : Young's modulus distribution for bar 1
%           -- (nDivisions+1:2*nDivisions) : Cross-sectional area distribution for bar 1
%           -- (2*nDivisions+1:3*nDivisions) : Young's modulus distribution for bar 2
%           -- (3*nDivisions+1:4*nDivisions) : Cross-sectional area distribution for bar 2
%       - SYSTEMPARAMETER : (1,3) double
%           -- (1) : length of bar 1
%           -- (2) : length of bar 2
%           -- (3) : number of interpolation points
%
%   Outputs:
%       - PERFORMANCEMEASURE : (nSample,1) double
%
%   See also interp1, trapz.
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

    length1 = systemParameter(1);
    length2 = systemParameter(2);
    nInterpolationPoints = systemParameter(3);

    interpolationMethod = 'pchip';

    nSample = size(designSample,1);
    nDivisions = size(designSample,2)/4;

    length1Base = linspace(0,length1,nDivisions);
    length1Fine = linspace(0,length1,nInterpolationPoints);

    length2Base = linspace(0,length2,nDivisions);
    length2Fine = linspace(0,length2,nInterpolationPoints);

    performanceMeasure = nan(nSample,1);
    for i=1:nSample
        youngsModulus1Base = designSample(i,1:nDivisions);
        area1Base = designSample(i,nDivisions+1:2*nDivisions);
        
        youngsModulus2Base = designSample(i,2*nDivisions+1:3*nDivisions);
        area2Base = designSample(i,3*nDivisions+1:4*nDivisions);
        
        if(length(length1Base) == 1)
            youngsModulus1Fine = repmat(youngsModulus1Base,1,nInterpolationPoints);
            area1Fine = repmat(area1Base,1,nInterpolationPoints);
        else
            youngsModulus1Fine = interp1(length1Base,youngsModulus1Base,length1Fine,interpolationMethod);
            area1Fine = interp1(length1Base,area1Base,length1Fine,interpolationMethod);
        end
        
        if(length(length2Base) == 1)
            youngsModulus2Fine = repmat(youngsModulus2Base,1,nInterpolationPoints);
            area2Fine = repmat(area2Base,1,nInterpolationPoints);
        else
            youngsModulus2Fine = interp1(length2Base,youngsModulus2Base,length2Fine,interpolationMethod);
            area2Fine = interp1(length2Base,area2Base,length2Fine,interpolationMethod);
        end
        
        stiffness1 = 1/trapz(length1Fine,1./(youngsModulus1Fine.*area1Fine));
        stiffness2 = 1/trapz(length2Fine,1./(youngsModulus2Fine.*area2Fine));

        performanceMeasure(i) = stiffness1 + stiffness2;
    end
end


