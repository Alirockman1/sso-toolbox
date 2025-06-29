function [boundingBoxStrict,boundingBoxRelaxed] = design_bounding_box(designSample,label)
%DESIGN_BOUNDING_BOX Find the minimal/maximal bounding box of labeled designs
%	DESIGN_BOUNDING_BOX returns bounding boxes that includes all designs 
%	labeled as 'true'. 
%	In its 'strict' output, the box has at least one 'true' design in all its
%	edges, being the smallest possible bounding box.
%	In its 'relaxed' output, the bounding box has at least one 'false' design
%	in all its edges, which is found by extending the 'strict' edges until a 
%	'false' design is found. 
%	Note, however, that that does not mean the 'relaxed' output will have no 
%	'false' designs inside itself, as 'false' designs in the corners between the
%	two boxes would not be detected by this expansion strategy. If higher 
%	purity of the relaxed bounding box is desired, a trimming strategy (with 
%	defined preferences) may become necessary.
%
%	BOUNDINGBOXSTRICT = DESIGN_BOUNDING_BOX(DESIGNSAMPLE) receives the
%	design sample points DESIGNSAMPLE and returns the strict bounding design 
%	box BOUNDINGBOXSTRICT, in other words, the smallest possible bounding box
%	that contains all those designs.
%
%	BOUNDINGBOXSTRICT = DESIGN_BOUNDING_BOX(DESIGNSAMPLE,LABEL) receives the
%	the boolean labels of each design in LABEL, where only designs labeled 
%	as 'true' get considered for the strict bounding design box. If no argument 
%   is given, it is assumed all designs must be in the bounding box (in other 
%	words, all are given the 'true' label).
%
%	[BOUNDINGBOXSTRICT,BOUNDINGBOXRELAXED] = DESIGN_BOUNDING_BOX(...) also
%	returns the relaxed bounding design box BOUNDINGBOXRELAXED, which is an 
%	expansion of the the strict box until a 'false'-labeled design is found in 
%	each edge.
%
%	Example
%		Create two random samples, one bad and one good, and visualize
%		the two bounding boxes.
%			designSampleGood = rand(100,2);
% 			designSampleBad = rand(100,2);
%
%			% make good designs in [-1,1] range, while bad designs in [-2,2] range
% 			designSampleGood = -1 + designSampleGood*2;
% 			designSampleBad = -2 + designSampleBad*4;
%
%			% concatenate arrays and label correctly
% 			designSample = [designSampleGood;designSampleBad];
% 			label = [true(size(designSampleGood,1),1);false(size(designSampleBad,1),1)];
%
%			% find bounding boxes
% 			[boundingBoxStrict,boundingBoxRelaxed] = DESIGN_BOUNDING_BOX(designSample,label);
%			
%			% visualize bounding boxes
%			figure;
% 			plot(designSampleGood(:,1),designSampleGood(:,2),'g.','MarkerSize',10);
% 			hold all;
% 			plot(designSampleBad(:,1),designSampleBad(:,2),'r.','MarkerSize',10);
% 			plot_design_box_2d(gcf,boundingBoxStrict,'EdgeColor','c');
% 			plot_design_box_2d(gcf,boundingBoxRelaxed,'EdgeColor','m');
% 			grid minor;
%		In the resulting figure, you can visually confirm the behavior of the
%		bounding boxes.
%
%	Inputs:
%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
%		- LABEL : (nSample,1) logical, optional
%	
%	Outputs:
%		- BOUNDINGBOXSTRICT : (2,nDesignVariable) double
%		- BOUNDINGBOXRELAXED : (2,nDesignVariable) double
%
%	See also box_apply_leanness, plot_design_box_2d.
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

    if(nargin<2)
        label = true(size(designSample,1),1);
    end

	nDimension = size(designSample,2);
	boundingBoxStrict = nan(2,nDimension);
	boundingBoxStrict(1,:) = min(designSample(label,:),[],1);
	boundingBoxStrict(2,:) = max(designSample(label,:),[],1);
	
	if(nargout>1)
		limitLowerBoundary = min(designSample,[],1);
		limitUpperBoundary = max(designSample,[],1);

		distanceLowerBound = boundingBoxStrict(1,:) - designSample;
		distanceUpperBound = designSample - boundingBoxStrict(2,:);
		belowLowerBound = (distanceLowerBound>0);
		aboveUpperBound = (distanceUpperBound>0);

		dimensionArray = 1:nDimension;

		for i=1:nDimension
			otherDimension = dimensionArray;
			otherDimension(i) = [];

			isInProjectedBox = is_in_design_box(designSample(:,otherDimension),boundingBoxStrict(:,otherDimension));

			lowerBoundaryCandidate = isInProjectedBox & belowLowerBound(:,i);
			newLowerBoundary = max(designSample(lowerBoundaryCandidate,i),[],1);
			boundingBoxRelaxed(1,i) = ternary(~isempty(newLowerBoundary),newLowerBoundary,limitLowerBoundary(i));

			upperBoundaryCandidate = isInProjectedBox & aboveUpperBound(:,i);
			newUpperBoundary = min(designSample(upperBoundaryCandidate,i),[],1);
			boundingBoxRelaxed(2,i) = ternary(~isempty(newUpperBoundary),newUpperBoundary,limitUpperBoundary(i));
		end
	end
end