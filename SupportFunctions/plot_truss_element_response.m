function figureHandle = plot_truss_element_response(nodePosition,nodeElement,elementResponse,varargin)
%PLOT_TRUSS_ELEMENT_RESPONSE Visualization of a response of each element
%   PLOT_TRUSS_ELEMENT_RESPONSE allows for the visualization of a given response
%   for each eleement of a truss, where the elements are colored according to
%   said response. This could be used to view axial forces, stresses, strains, 
%   and similar quantities.
%
%   FIGUREHANDLE = PLOT_TRUSS_ELEMENT_RESPONSE(NODEPOSITION,NODEELEMENT,
%   ELEMENTRESPONSE) plots a truss defined by its nodal positions NODEPOSITION,
%   element definitions NODEELEMENT and given response ELEMENTRESPONSE in a new
%   figure FIGUREHANDLE. Each element is colored according to its respective
%   entry in ELEMENTRESPONSE.
%
%   FIGUREHANDLE = PLOT_TRUSS_ELEMENT_RESPONSE(NODEPOSITION,NODEELEMENT,
%   ELEMENTRESPONSE,ELEMENTCROSSSECTIONAREA) uses the cross-section areas of 
%   each truss and plots each element with smaller / larger linewidth 
%   accordingly. By default if this is not given, the maximum linewidth is used.
%
%   FIGUREHANDLE = PLOT_TRUSS_ELEMENT_RESPONSE(...,NAME,VALUE,...) allows for  
%   the specification of additional options. These are:
%       - 'MinimumLinewidth' : if empty, scaling of elements is done assuming
%       the conversion is centered around the origin (as in, an element of zero 
%       area would have zero linewidth). If not empty, this is the linewidth 
%       used in the plot for the truss elements with the smallest area. Default:
%       empty.
%       - 'MaximumLinewidth' : linewidth used in the plot for the truss elements
%       with the largest area (and also for trusses where all areas are equal)
%       Default: 1.
%       - 'PlotOptions' : additional options to be used for 'plot'/'plot3' when
%       plotting each element of the truss (both undeformed and deformed). 
%       'Color' and 'linewidth' should not be set here, but instead, in the 
%       previous options. By default, 'Marker' is set to 'o', but that can be
%       overwritten, and other options can be added. 
%
%   Input:
%       - NODEPOSITION : (nNode,2) double OR (nNode,3) double 
%       - NODEELEMENT : (nElement,2) integer 
%       - ELEMENTCROSSSECTIONAREA : double OR (nElement,1) double
%       - 'MinimumLinewidth' : double
%       - 'MaximumLinewidth' : double
%       - 'PlotOptions' : (1,nOptions) cell
%
%   Output:
%       - FIGUREHANDLE : figure_handle
%
%   See also truss_analysis, plot_truss_deformation.
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

	parser = inputParser;
    parser.addOptional('ElementCrossSectionArea',1);
    parser.addParameter('MinimumLinewidth',[],@(x)isnumeric(x));
    parser.addParameter('MaximumLinewidth',1,@(x)isnumeric(x));
    parser.addParameter('PlotOptions',{},@(x)iscell(x));
	parser.parse(varargin{:});
	options = parser.Results;

    nElement = size(nodeElement,1);
    elementCrossSectionArea = parser.Results.ElementCrossSectionArea;
    if(length(elementCrossSectionArea)==1)
        elementCrossSectionArea = elementCrossSectionArea * ones(nElement,1);
    end

	defaultPlotOptions = {'Marker','o'};
    [~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,options.PlotOptions);

    % calculate linewidth based on areas
    if(isempty(options.MinimumLinewidth) || isnan(options.MinimumLinewidth))
        minArea = 0;
        minLinewdith = 0;
    else
        minArea = min(elementCrossSectionArea);
        minLinewdith = options.MinimumLinewidth;
    end

    % linearly interpolate between minimum and maximum values
    maxArea = max(elementCrossSectionArea);
    rangeArea = maxArea - minArea;

    maxLinewidth = options.MaximumLinewidth;
    rangeLinewidth = maxLinewidth - minLinewdith;
    
    elementLinewidth = minLinewdith + rangeLinewidth.*(elementCrossSectionArea-minArea)./rangeArea;

    % begin plotting truss elements
    nElement = size(nodeElement,1);
	nDimension = size(nodePosition,2);

	figureHandle = figure;
	hold all;
	colorbar;
    clim([min(elementResponse),max(elementResponse)]);
    colorSteps = colormap;
    colorIndex = round(interp1([min(elementResponse);max(elementResponse)],[1 size(colorSteps,1)],elementResponse));
	for i=1:nElement
        trussX = nodePosition(nodeElement(i,:),1);
    	trussY = nodePosition(nodeElement(i,:),2);

        if(nDimension==2)
            plotInputArgument = {trussX,trussY};
            plotFunction = @plot;
        else
            trussZ = nodePosition(nodeElement(i,:),3);
            plotInputArgument = {trussX,trussY,trussZ};
            plotFunction = @plot3;
        end

        plotFunction(plotInputArgument{:},...
            'Color',colorSteps(colorIndex(i),:),...
            'Linewidth',elementLinewidth(i),...
            plotOptions{:});
    end
end