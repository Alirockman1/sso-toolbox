function activate_graphics_object(graphicsHandle)
%ACTIVATE_GRAPHICS_OBJECT Activates a MATLAB graphics object based on its type
%   ACTIVATE_GRAPHICS_OBJECT makes a graphics object the current active object.
%   It works with different types of MATLAB graphics objects such as figures,
%   axes, UI figures, etc.
%
%   Input:
%       - GRAPHICSHANDLE : handle
%
%   Examples:
%       fig = figure();
%       activate_figure_axis(fig);
%
%       ax = axes();
%       activate_figure_axis(ax);
%
%       uifig = uifigure();
%       activate_figure_axis(uifig);
%
%   See also figure, axes, uifigure.
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

    % Input validation
    if ~ishandle(graphicsHandle) && ~isgraphics(graphicsHandle)
        error('Input must be a valid graphics handle');
    end
    if numel(graphicsHandle) ~= 1
        error('Input must be a scalar graphics handle.');
    end

    % Get the type of the graphics object
    objectType = lower(get(graphicsHandle, 'Type'));

    % Activate the parent figure, if available.
    parentFig = ancestor(graphicsHandle, 'figure');
    if ~isempty(parentFig)
        try
            % For legacy figures, this brings the figure to the front.
            figure(parentFig);
        catch
            % In case of errors (e.g., with uifigure), ensure visibility.
            if isprop(parentFig, 'Visible')
                parentFig.Visible = 'on';
            end
            drawnow;
        end
    end

    % Activate the object based on its type
    switch objectType
        case {'figure', 'uifigure'}
            % For figure objects, attempt to bring it forward.
            try
                figure(graphicsHandle);
            catch
                % For uifigure, if figure() is not supported, force visibility.
                if isprop(graphicsHandle, 'Visible')
                    graphicsHandle.Visible = 'on';
                end
            end
            
        case 'axes'
            % Set the current axes.
            axes(graphicsHandle);
            
        case 'uiaxes'
            % For UI axes, activate by setting the 'Selected' property if it exists.
            if isprop(graphicsHandle, 'Selected')
                graphicsHandle.Selected = 'on';
            else
                warning('The provided UI axes does not support activation via the ''Selected'' property.');
            end
            
        case 'uitab'
            % For a UI tab, activate it using the 'Selected' property.
            if isprop(graphicsHandle, 'Selected')
                graphicsHandle.Selected = true;
            else
                warning('The provided UI tab does not support activation via the ''Selected'' property.');
            end
            
        otherwise
            % For other graphics object types, issue a warning.
            warning('Unable to activate the specified graphics object type: %s', objectType);
    end
end
