function trimmingOrder = trimming_order(isExclude,score,varargin)
%TRIMMING_ORDER Define the order of designs to remove in the trimming operation
%	TRIMMING_ORDER allows for the determination of the order that designs will
%	be removed during the trimming operation. This is determined either randomly
%	or with the respective score of each design that is not acceptable.
%
%	TRIMMINGORDER = TRIMMING_ORDER(ISEXCLUDE,SCORE) will order the designs
%	labeled as not acceptable ('false') in ISEXCLUDE from lowest-to-highest 
%	SCORE, returning that order in TRIMMINGORDER.
%
%	TRIMMINGORDER = TRIMMING_ORDER(...NAME,VALUE,...) allows for the choice
%	of order to use. The option name 'Order' can assume the following values:
%		- 'score-low-to-high' / 'lth' / 'low-to-high' / 'lowest-first': order 
%		them from lowest score first, therefore trimming first closest to the 
%		boundary.
%		- 'score-high-to-low' / 'htl' / 'high-to-low' / 'highest-first' : order 
%		them from highest score first, therefore trimming first furthest to the 
%		boundary.
%		- 'score-both' / 'score' : both orders for the score, going both from
%		low to high scores, and after that high to low scores in a second pass.
%		- 'random' / 'rand' : randomly pick the order of elimination.
%		- a function handle in the form 'orderExclude = f(scoreExclude)', which
%		can return an order given all the scores of designs that must be 
%		excluded.
%
%   Inputs:
%       - ISEXCLUDE : (nSample,1) logical
%       - SCORE : (nSample,1) double
%		- 'Order' : char OR string OR function_handle 
%   
%   Output:
%       - TRIMMINGORDER : (nExclude,nOrder) integer
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

	tagScoreLowToHigh = {'score-low-to-high','lth','low-to-high','lowest-first'};
	tagScoreHighToLow = {'score-high-to-low','htl','high-to-low','highest-first'};
	tagScoreBoth = {'score-both','score'};
	tagRandom = {'random','rand'};

	parser = inputParser;
	parser.addRequired('isExclude',@(x)islogical(x)&&size(x,2)==1);
	parser.addRequired('DesignScore',@(x)isnumeric(x)&&size(x,2)==1);
	parser.addParameter('OrderPreference','score',@(x)any(strcmpi(x,[tagScoreLowToHigh,tagScoreHighToLow,tagRandom,tagScoreBoth]))||isa(x,'function_handle'));
	parser.parse(isExclude,score,varargin{:});

	orderPreference = parser.Results.OrderPreference;

	% find designs that must be excluded
	scoreExclude = score(isExclude);
	nExclude = size(scoreExclude,1);

    % if nothing is to be trimmed, return
    if(nExclude<=0)
        trimmingOrder = [];
        return;
    end
	
	% define order of trimming
	if(isa(orderPreference,'function_handle'))
		orderExclude = orderPreference(scoreExclude);
    elseif(any(strcmpi(orderPreference,tagScoreLowToHigh))) 
        [~,orderExclude] = sort(scoreExclude,'ascend');
    elseif(any(strcmpi(orderPreference,tagScoreHighToLow)))
        [~,orderExclude] = sort(scoreExclude,'descend');
    elseif(any(strcmpi(orderPreference,tagScoreBoth)))
    	orderExclude = nan(nExclude,2);
        [~,orderExclude(:,1)] = sort(scoreExclude,'ascend');
        orderExclude(:,2) = flip(orderExclude(:,1));
    elseif(strcmpi(orderPreference,'random'))
        orderExclude = randperm(nExclude);
    else
    	orderExclude = (1:nExclude)';
    end

    % convert to all-elements base index
    trimmingOrder = nan(size(orderExclude));
    for i=1:size(orderExclude,2)
        trimmingOrder(:,i) = convert_index_base(isExclude,orderExclude(:,i),'backward');
    end
end