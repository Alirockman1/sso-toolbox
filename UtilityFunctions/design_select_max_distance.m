function designSelection = design_select_max_distance(referenceDesign,queryDesign,nSelect,varargin)
%DESIGN_SELECT_MAX_DISTANCE Select designs with maximum distance to reference
%	DESIGN_SELECT_MAX_DISTANCE returns an ordered index of the designs in a 
%	query that are the most distant from reference designs and from previous
%	selections on that query. 
%
%	DESIGNSELECTION = DESIGN_SELECT_MAX_DISTANCE(REFERENCEDESIGN,QUERYDESIGN)
%	chooses from the query QUERYDESIGN the designs most distant from 
%	REFERENCEDESIGN and from previously selected designs, returning their 
%	indices in DESIGNSELECTION. This means the design in the first index of
%	DESIGNSELECTION is the most distant from all in REFERENCEDESIGN; the second
%	is the most distant from all designs REFERENCEDESIGN and the first design
%	selected from QUERYDESIGN; and so on. In this case, default options are used
%	for 'knnsearch', and it is assumed all query designs should be ordered in
%	this way. If REFERENCEDESIGN is empty, the first design of QUERYDESIGN
%	is used to start the process.
%
%	DESIGNSELECTION = DESIGN_SELECT_MAX_DISTANCE(REFERENCEDESIGN,QUERYDESIGN,
%	NSELECT) allows one to set how many designs NSELECT should be ordered in  
%	this manner. By default, they are  all ordered.
%
%	DESIGNSELECTION = DESIGN_SELECT_MAX_DISTANCE(REFERENCEDESIGN,QUERYDESIGN,
%	NSELECT,...NAME,VALUE,...) also allows for setting options for the distance 
%	computation through 'knnsearch' with the name-value pair arguments.
%
%   Inputs:
%       - REFERENCEDESIGN : (nReference,nDesignVariable) double 
%       - QUERYDESIGN : (nQuery,nDesignVariable) double
%       - NSELECT :  integer
%		- Name-value pair arguments : passed directly to 'knnsearch'.
%
%   Outputs:
%       - DESIGNSELECTION : (nSelect,1) double
%
%   See also knnsearch, design_sort_min_distance.
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

	if(nargin<3 || isempty(nSelect))
		nSelect = size(queryDesign,1);
	end

	% initialize 
	nDesignVariable = size(queryDesign,2);
	designSelection = nan(nSelect,1);

	if(isempty(referenceDesign))
		% if no reference is given, start from the first design in the query
		% and start calculating the distances from there (max(inf array) returns
		% the first index)
		distance = inf(size(queryDesign,1),1);
	else
		% start by finding distances between the references and query
		[~,distance] = knnsearch(referenceDesign,queryDesign,varargin{:});
	end

	% select the design with greatest distance and update the distance array
	% with the distance relative to the chosen design
	for i=1:nSelect
		% get the queryDesign with the greatest distance to obtained samples
		[~,designSelection(i)] = max(distance);

		% update distances with new entry, assuming it will be evaluated
		[~,newDistance] = knnsearch(queryDesign(designSelection(i),:),queryDesign,varargin{:});
		distance = min(distance,newDistance);
	end
end