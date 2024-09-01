classdef (Abstract) BottomUpMappingBase
%BOTTOMUPMAPPINGBASE Bottom-Up Mapping abstract class
%	BOTTOMUPMAPPINGBASE provides the necessary interfaces for particular 
%   implementations that may serve as bottom-up mappings in the context of
%	systems design. The goal of this class and its derivations is to receive 
%   as input design sample points and return the system response.
%   Additionally, if physical feasibility is being considered, the indicators or
%   measures of said physical feasibility should also be returned.
%
% 	As an abstract class, particular implementations have to be created using 
%   this as base for one to be able to create objects.
%
%	BOTTOMUPMAPPINGBASE methods:
%		- response : a function that receives the design sample points and  
%       returns the system response for those variables in terms of performance
%       and physical feasibility.
%
%   See also BottomUpMappingFunction, BottomUpMappingPython.
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

	methods(Abstract) 
		%RESPONSE Getting system response based on the design sample points
		%	RESPONSE is the main function of any bottom-up mapping, where design sample 
		%	points are given and the system response in terms of performance and  
		%	physical feasibility is returned.
        %
        %   PERFORMANCEMEASURE = OBJ.RESPONSE(DESIGNSAMPLE) returns the performance
        % 	measure of each given design in PERFORMANCEMEASURE. This is returned as an 
        %	array with each row being correspondent to the respective design.
        %
        %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = OBJ.RESPONSE(...) also 
        %	returns the physical feasibility measure of each given design in 
        %	PHYSICALFEASIBILITYMEASURE; similar to PERFORMANCEMEASURE, this is given as 
        %	an array with each row being correspondent to the respective design.
        %
        %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE,SYSTEMOUTPUT] = 
        %	OBJ.RESPONSE(...) also returns any more information from the system response 
        %	computation in SYSTEMOUTPUT; what that is and its type is entirely dependent 
        %   on the specific class implementation.
        %
        %   Inputs:
        %       - OBJ : BottomUpMappingBase
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - PERFORMANCEMEASURE : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - SYSTEMOUTPUT : class-defined
        %
        %   See also BottomUpMappingBase.
		[performanceMeasure,physicalFeasibilityMeasure,systemOutput] = response(obj,designSample)
	end
end