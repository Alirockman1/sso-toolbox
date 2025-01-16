function [handleUndeformed,handleDeformed] = plot_truss_deformation(figureHandle,nodePosition,nodeElement,varargin)
%PLOT_TRUSS_DEFORMATION Visualization of a Truss (and its deformed state)
%   PLOT_TRUSS_DEFORMATION plots a truss structure in the given figure, and if
%   the displacement of the nodes is given, also its deformed state. 
%
%   HANDLEUNDEFORMED = PLOT_TRUSS_DEFORMATION(FIGUREHANDLE,NODEPOSITION,
%   NODEELEMENT) plots a truss defined by its nodal positions NODEPOSITION and
%   element definitions NODEELEMENT onto the figure FIGUREHANDLE, returning the
%   handle for the line object HANDLEUNDEFORMED.
%
%   [HANDLEUNDEFORMED,HANDLEDEFORMED] = PLOT_TRUSS_DEFORMATION(FIGUREHANDLE,
%   NODEPOSITION,NODEELEMENT,NODEDISPLACEMENT) also plots the deformed truss
%   based on the nodal displacements NODEDISPLACEMENT, and returns the handle
%   for that line object in HANDLEDEFORMED.
%
%   [HANDLEUNDEFORMED,HANDLEDEFORMED] = PLOT_TRUSS_DEFORMATION(FIGUREHANDLE,
%   NODEPOSITION,NODEELEMENT,NODEDISPLACEMENT,ELEMENTCROSSSECTIONAREA) uses
%   the cross-section areas of each truss and plots each element with smaller / 
%   larger linewidth accordingly. By default if this is not given, the maximum 
%   linewidth is used.
%
%   [...] = PLOT_TRUSS_DEFORMATION(...,NAME,VALUE,...) allows for the 
%   specification of additional options. These are:
%       - 'ColorUndeformed' : color of the truss when undeformed. Default: 'b'.
%       - 'ColorDeformed' : color of the truss when deformed. Default: 'r'.
%       - 'DisplacementScaleFactor' : a scale factor for the displacements when
%       plotting the deformed truss, which may allow for better visualization
%       of it if the displacement is small. Default: 1.
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
%       - FIGUREHANDLE : figure_handle
%       - NODEPOSITION : (nNode,2) double OR (nNode,3) double 
%       - NODEELEMENT : (nElement,2) integer
%       - NODEDISPLACEMENT : (nNode,2) double OR (nNode,3) double 
%       - ELEMENTCROSSSECTIONAREA : double OR (nElement,1) double
%       - 'ColorUndeformed' : color (char OR string OR (1,3) double)
%       - 'ColorDeformed' : color (char OR string OR (1,3) double)
%       - 'DisplacementScaleFactor' : double
%       - 'MinimumLinewidth' : double
%       - 'MaximumLinewidth' : double
%       - 'PlotOptions' : (1,nOptions) cell
%
%   Output:
%       - HANDLEUNDEFORMED : line object
%       - HANDLEDEFORMED : line object
%
%   See also truss_analysis, plot_truss_element_response.
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

	parser = inputParser;
    parser.addOptional('NodeDisplacement',[]);
    parser.addOptional('ElementCrossSectionArea',1);
    parser.addParameter('DisplacementScaleFactor',1);
    parser.addParameter('MinimumLinewidth',[],@(x)isnumeric(x));
    parser.addParameter('MaximumLinewidth',1,@(x)isnumeric(x));
    parser.addParameter('TrussPlotOptions',{},@(x)iscell(x));
    parser.addParameter('UndeformedPlotOptions',{});
    parser.addParameter('DeformedPlotOptions',{});
	parser.parse(varargin{:});
    options = parser.Results;

    nodeDisplacement = parser.Results.NodeDisplacement;

    nElement = size(nodeElement,1);
    elementCrossSectionArea = parser.Results.ElementCrossSectionArea;
    if(length(elementCrossSectionArea)==1)
		elementCrossSectionArea = elementCrossSectionArea * ones(nElement,1);
	end

    defaultTrussPlotOptions = {'Color',color_palette_tol('blue'),'Marker','o'};
    [~,trussPlotOptions] = merge_name_value_pair_argument(defaultTrussPlotOptions,options.TrussPlotOptions);

    defaultUndeformedPlotOptions = {'Linestyle','-'};
    [~,undeformedPlotOptions] = merge_name_value_pair_argument(defaultUndeformedPlotOptions,trussPlotOptions,options.UndeformedPlotOptions);
    
    defaultDeformedPlotOptions = {'Linestyle',':'};
    [~,deformedPlotOptions] = merge_name_value_pair_argument(defaultDeformedPlotOptions,trussPlotOptions,options.DeformedPlotOptions);
    
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
    nDimension = size(nodePosition,2);
    handleDeformed = [];
    isPlotDeformed = (~isempty(nodeDisplacement));

    figure(figureHandle);
    hold all;
    for i=1:nElement
    	undeformedTrussX = nodePosition(nodeElement(i,:),1);
    	undeformedTrussY = nodePosition(nodeElement(i,:),2);
        if(isPlotDeformed)
    	   deformedTrussX = nodePosition(nodeElement(i,:),1) + options.DisplacementScaleFactor*nodeDisplacement(nodeElement(i,:),1);
    	   deformedTrussY = nodePosition(nodeElement(i,:),2) + options.DisplacementScaleFactor*nodeDisplacement(nodeElement(i,:),2);
        end

        if(nDimension==2)
            plotInputArgumentUndeformed = {undeformedTrussX,undeformedTrussY};

            if(isPlotDeformed)
                plotInputArgumentDeformed = {deformedTrussX,deformedTrussY};
            end

            plotFunction = @plot;
        else
            undeformedTrussZ = nodePosition(nodeElement(i,:),3);
            plotInputArgumentUndeformed = {undeformedTrussX,undeformedTrussY,undeformedTrussZ};

            if(isPlotDeformed)
                deformedTrussZ = nodePosition(nodeElement(i,:),3) + options.DisplacementScaleFactor*nodeDisplacement(nodeElement(i,:),3);
                plotInputArgumentDeformed = {deformedTrussX,deformedTrussY,deformedTrussZ};
            end

            plotFunction = @plot3;
        end

        handleUndeformed = plotFunction(plotInputArgumentUndeformed{:},...
            'Linewidth',elementLinewidth(i),...
            undeformedPlotOptions{:});
        if(isPlotDeformed)
            handleDeformed = plotFunction(plotInputArgumentDeformed{:},...
                'Linewidth',elementLinewidth(i),...
                deformedPlotOptions{:});
        end
    end
end