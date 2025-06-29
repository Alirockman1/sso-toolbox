%
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

%% Cleanup
close all;
fclose all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(18);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% define ellipsoid - 2d
ellipsoidCenter2d = [2 1];
ellipsoidMajorSemiAxis2d = [4 3];
ellipsoidRotationAngle2d = 0.5;


%% design space to explore
designSpace2d = [-3 -3 ; 7 5];
nTest = 1000;


%%
designSample = sampling_latin_hypercube(designSpace2d,1000);
performanceMeasure2d = ellipse_2d(designSample,[ellipsoidCenter2d,ellipsoidMajorSemiAxis2d,ellipsoidRotationAngle2d]);
label2d = (performanceMeasure2d<=1);

figure;
hold all;
plot_ellipse_2d(gcf,ellipsoidCenter2d,ellipsoidMajorSemiAxis2d,ellipsoidRotationAngle2d,...
	'PatchOptions',{'FaceAlpha',0});
plot(designSample(label2d,1),designSample(label2d,2),'g.');
plot(designSample(~label2d,1),designSample(~label2d,2),'r.');
xlabel('x');
ylabel('y');
axis('equal');
grid minor;


%% define ellipsoid - 3d
ellipsoidCenter3d = [2 1 3];
ellipsoidMajorSemiAxis3d = [4 3 5];
ellipsoidRotationAngle3d = [0.5 -0.1 0.2];


%% design space to explore
designSpace3d = [-3 -3 -3; 7 5 9];
nTest = 1000;


%%
designSample = sampling_latin_hypercube(designSpace3d,1000);
performanceMeasure3d = ellipsoid_3d(designSample,[ellipsoidCenter3d,ellipsoidMajorSemiAxis3d,ellipsoidRotationAngle3d]);
label3d = (performanceMeasure3d<=1);

figure;
hold all;
plot_ellipsoid_3d(gcf,ellipsoidCenter3d,ellipsoidMajorSemiAxis3d,ellipsoidRotationAngle3d,...
	'SurfOptions',{'FaceAlpha',0});
plot3(designSample(label3d,1),designSample(label3d,2),designSample(label3d,3),'g.');
plot3(designSample(~label3d,1),designSample(~label3d,2),designSample(~label3d,3),'r.');
xlabel('x');
ylabel('y');
zlabel('z');
axis('equal');
grid minor;


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

