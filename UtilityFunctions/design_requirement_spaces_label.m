function [isAcceptable,isUseful] = design_requirement_spaces_label(requirementSpacesType,isGoodPerformance,isPhysicallyFeasible)
%DESIGN_REQUIREMENT_SPACES_LABEL Acceptable/useful labels for requirement spaces
%   DESIGN_REQUIREMENT_SPACES_LABEL transforms labels in terms of whether or not  
%   performance requirements are met and whether or not a design is physically 
%   feasible into labels that indicate if each given design is acceptable or not 
%   and useful or not, according to the definitions of the chosen problem type
%   of requirement spaces. The definitions are:
%       - 'Omega0' : acceptable designs are good performance, all designs are
%       useful.
%       - 'Omega1' : acceptable designs are good performance and physically 
%       feasible, all designs are useful.
%       - 'Omega2' : acceptable designs are good performance, useful designs
%        are physically feasible.
%       - 'Omega3' : acceptable designs are good performance or physically 
%       infeasible, useful designs are physically feasible.
%
%   ISACCEPTABLE = DESIGN_REQUIREMENT_SPACES_LABEL(REQUIREMENTSPACESTYPE,
%   ISGOODPERFORMANCE,ISPHYSICALLYFEASIBLE) can, for requirement spaces type
%   REQUIREMENTSPACESTYPE, use the labels of each design (regarding whether it
%   meets the performance requirements ISGOODPERFORMANCE and whether it is 
%   physically feasible or not ISPHYSICALLYFEASIBLE) to determine if each 
%   design is acceptable or not ISACCEPTABLE. By default, the type of 
%   requirement spaces that is assumed is 'Omega2'. As for the labels, 'true'  
%   indicates that design sample meets all performance requirements / is 
%   physically feasible / is acceptable, and 'false' indicates the opposite.
%
%   [ISACCEPTABLE,ISUSEFUL] = DESIGN_REQUIREMENT_SPACES_LABEL(...) also returns
%   a label for indicating whether each design is useful or not ISUSEFUL; like
%   before, 'true' indicates that design is useful.
%
%   Inputs:
%       - REQUIREMENTSPACESTYPE : string 
%       - ISGOODPERFORMANCE : (nSample,nTest) logical
%       - ISPHYSICALLYFEASIBLE :  (nSample,nTest) logical
%
%   Outputs:
%       - ISACCEPTABLE : (nSample,nTest) logical
%       - ISUSEFUL : (nSample,nTest) logical
%
%   See also sso_box_stochastic, sso_component_stochastic.
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

    if(strcmpi(requirementSpacesType,'Omega0'))
        isAcceptable = isGoodPerformance; 
        isUseful = true(size(isPhysicallyFeasible)); % design space
    elseif(strcmpi(requirementSpacesType,'Omega1'))
        isAcceptable = isGoodPerformance & isPhysicallyFeasible; 
        isUseful = true(size(isPhysicallyFeasible)); % design space
    elseif(strcmpi(requirementSpacesType,'Omega3'))
        isAcceptable = isGoodPerformance | (~isPhysicallyFeasible);
        isUseful = isPhysicallyFeasible; 
    else%if(strcmpi(requirementSpacesType,'Omega2'))
        isAcceptable = isGoodPerformance;
        isUseful = isPhysicallyFeasible; 
    end
end