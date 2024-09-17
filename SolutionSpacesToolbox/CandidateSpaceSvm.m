classdef CandidateSpaceSvm < CandidateSpaceBase
%CandidateSpaceSvm Candidate Space defined by a Support Vector Machine
%   CandidateSpaceSvm can be used to define candidate spaces for  
%   component solution space computation. In this case, the candidate space is 
%   defined by a SVM trained on the labeled data given for design samples.
%
%   This class is derived from CandidateSpace.
%
%   CandidateSpaceSvm properties:
%       - DesignSampleDefinition : design samples used in the candidate space
%       definition.
%       - LabelDefinition : labels of the design samples used in the candidate
%       space definition.
%       - SamplingBox : bounding box around the positive region of the candidate
%       space that can be used to better attempt to sample inside said region.
%       - Measure : measure of the candidate space.
%       - ActiveDesigns : design samples used in definition with positive label.
%       - SamplingBoxSlack : factor to determine how strict/relaxed the sampling
%       box is.
%       - Svm : support vector machine used to define candidate space.
%       - SvmBalanceCost : option on how to weight false positives in training.
%       - SvmTrainingOptions : options for the training of the SVM.
%       - GrowthDistanceOptions : options for the calculation of distance during
%       the candidate space growth operation.
%
%   CandidateSpaceSvm methods:
%       - define_candidate_space : create a candidate space based on design 
%       samples that are labeled as inside/outside.
%       - grow_candidate_space : expand the candidate space by a given factor.
%       - is_in_candidate_space : verify if given design samples are inside 
%       the candidate space.
%
%   See also ClassificationSVM.
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

    properties(SetAccess = protected, Dependent)
        %DESIGNSAMPLEDEFINITION Design sample points used in candidate space definition
        %   DESIGNSAMPLEDEFINITION are the sample points used in the definition of the 
        %   current candidate space.
        %
        %   DESIGNSAMPLEDEFINITION : (nSample,nDesignVariable) double
        %
        %   See also IsInsideDefinition, IsShapeDefinition.
        DesignSampleDefinition
        
        %ISINSIDEDEFINITION Labels of sample points used in candidate space definition
        %   ISINSIDEDEFINITION are the labels of design samples used in the definition 
        %   of the current candidate space. A label of 'true' indicates the respective 
        %   design point is inside the candidate space.
        %
        %   ISINSIDEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsShapeDefinition.
        IsInsideDefinition

        %ISSHAPEDEFINITION Labels if sample points from definition contributes to shape
        %   ISSHAPEDEFINITION is a logical array where 'true' values indicate that that
        %   design point (in the respective row) from the definition sample actively 
        %   contributes to the shape of the candidate space. For example, for a
        %   convex hull, this would be the convex hull index points.
        %
        %   ISSHAPEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        IsShapeDefinition

        %MEASURE Size measure of the candidate space
        %   MEASURE is a value that works as the measure of the candidate space.  
        %   For CandidateSpaceSvm, MEASURE is the volume of the inside region which is
        %   estimated with the Monte Carlo integration method.
        %
        %   MEASURE : double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition. 
        Measure
    end

    properties (SetAccess = protected)
        %SVM Support Vector Machine model
        %   SVM is the object that contains the support vector machine used for
        %   the definition and classification of the candidate space.
        %
        %   SVM : ClassificationSVM 
        %
        %   See also fitcsvm, ClassificationSVM.
        Svm

        %SVMTRAININGOPTIONS Options for SVM training during definition
        %   SVMTRAININGOPTIONS is a cell that contains further options for the
        %   training of the support vector machine. In its default value, it
        %   contains the following parameters:
        %       - Standardize : true
        %       - KernelFunction : 'gaussian'
        %       - Verbose : 0
        %   These and other options of 'fitcsvm' can be chosen as desired, with
        %   the following two exceptions:
        %       - 'Cost' should not be changed, as it is already defined (see 
        %       'InsidePriorityRatio' property).
        %       - 'ClassNames', which is fixed to [false,true].
        %   SVMTRAININGOPTIONS should be defined as a cell with name-value pair 
        %   arguments or a structure.
        %
        %   SVMTRAININGOPTIONS : cell with name-value pairs OR structure
        %
        %   See also fitcsvm, InsidePriorityRatio.
        SvmTrainingOptions

        %INSIDEPRIORITYRATIO Cost of false negative classification in SVM training
        %   INSIDEPRIORITYRATIO determines the weight of making false negative errors
        %   (is inside, classified outside) when training the Support Vector Machine 
        %   that defines the candidate space. The weight of making false positive errors
        %   (is outside, classified inside) is fixed at 1.
        %   INSIDEPRIORITYRATIO can be determined as either a fixed number or as a 
        %   function which takes the number of inside/outside points in the training 
        %   dataset.
        %   Default value is 1. In practice, with larger values, the candidate space
        %   becomes larger than it should (in other words, it includes more designs that
        %   should be outside), and with smaller values, it becomes smaller than it 
        %   should (designs that should be inside are judged as outside).
        %
        %   INSIDEPRIORITYRATIO : double OR function_handle with the format
        %   'insidePriorityRatio = f(nInsideTraining,nOutsideTraining)' 
        %
        %   See also fitcsvm, SvmTrainingOptions.
        InsidePriorityRatio

        %GROWTHDISTANCEOPTIONS Options when computing distances in the growth algorithm
        %   GROWTHDISTANCEOPTIONS are the extra options to be used when finding the
        %   distance between sample points in the growth algorithm. 
        %   By default in its base implementation, these would be options for the 
        %   'knnsearch' funciton.
        %
        %   GROWTHDISTANCEOPTIONS : (1,nOption) cell
        %
        %   See also grow_candidate_space, knnsearch.
        GrowthDistanceOptions
    end
    
    methods
        function obj = CandidateSpaceSvm(designSpaceLowerBound,designSpaceUpperBound,varargin)
        %CANDIDATESPACESVM Constructor
        %   CANDIDATESPACESVM is a constructor initializes an object of this class.
        %
        %   OBJ = CANDIDATESPACESVM(DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND) creates
        %   an object of the CandidateSpaceSvm class and sets its design space 
        %   boundaries to its input values DESIGNSPACELOWERBOUND and 
        %   DESIGNSPACEUPPERBOUND. Other properties are set to empty.
        %
        %   OBJ = CandidateSpaceDecisionTree(...,NAME,VALUE,...) also allows one to set
        %   specific options for the object. This can be 
        %       - 'SamplingBoxSlack' : where the boundaries of the sampling box will be
        %       relative to the strictest bounding box and the most relaxed bounding 
        %       box. A value of 0 means no slack and therefore the sampling box will be
        %       the most strict one possible, and 1 means the sampling box will be the 
        %       most relaxed.
        %       - 'SvmTrainingOptions' : options for the training of the support vector
        %       machine; see more on the respective property header.
        %       - 'InsidePriorityRatio' : weight factor definition for false negative 
        %       errors comitted during training. See more on the respective property
        %       header.
        %
        %   Inputs:
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %       - 'SamplingBoxSlack' : double
        %       - 'SvmTrainingOptions' : cell with name-value pairs OR structure
        %       - 'InsidePriorityRatio' : function handle OR double
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceSvm
        %   
        %   See also SvmTrainingOptions, InsidePriorityRatio.

            % parse inputs
            parser = inputParser;
            parser.addRequired('designSpaceLowerBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addRequired('designSpaceUpperBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addParameter('SamplingBoxSlack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
            parser.addParameter('SvmTrainingOptions',{},@(x)iscell(x));
            parser.addParameter('InsidePriorityRatio',1,@(x)isa(x,'function_handle')||isnumeric(x));
            parser.addParameter('GrowthDistanceOptions',{},@(x)iscell(x));
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            % initialize properties
            obj.Svm = [];
            obj.DesignSpaceLowerBound = parser.Results.designSpaceLowerBound;
            obj.DesignSpaceUpperBound = parser.Results.designSpaceUpperBound;
            obj.InsidePriorityRatio = parser.Results.InsidePriorityRatio;
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;
            obj.SvmTrainingOptions = parser.Results.SvmTrainingOptions;

            defaultGrowthDistanceOptions = {};
            [~,obj.GrowthDistanceOptions] = merge_name_value_pair_argument(...
                defaultGrowthDistanceOptions,parser.Results.GrowthDistanceOptions);
        end
        
        function obj = define_candidate_space(obj,designSample,isInside)
        %define_candidate_space Initial definition of the candidate space
        %   define_candidate_space uses labeled design samples to create the
        %   positive/negative regions of the candidate space. For 
        %   CandidateSpaceSvm, this means the labeled design samples are used
        %   to train a Support Vector Machine.
        %
        %   OBJ = OBJ.define_candidate_space(DESIGNSAMPLE,LABEL) receives the
        %   design samples in DESIGNSAMPLE and their positive/negative 
        %   (true/false) labels in LABEL, and returns a candidate space object  
        %   OBJ with the new definition.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceSvm
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - LABEL : (nSample,1) logical
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceSvm object
        %   
        %   See also fitcsvm, is_in_candidate_space.

            if(nargin<3)
                isInside = true(size(designSample,1),1);
            end

            defaultSvmTrainingOptions = {...
                'Standardize',true,...
                'KernelFunction','rbf',...
                'Verbose',0};

            % Svm Training
            if(all(isInside == isInside(1))) % single class
                obj.Svm = fitcsvm(designSample,isInside,defaultSvmTrainingOptions{:}); 
            else % two classes
                % Weighting factor of Type II error 
                % false negative - is positive (inside) design, predicted negative (outside) design
                if(isa(obj.InsidePriorityRatio,'function_handle'))
                    nOutsideTraining = sum(~isInside);
                    nInsideTraining = sum(isInside);

                    falseNegativeWeight = obj.InsidePriorityRatio(nInsideTraining,nOutsideTraining);
                else
                    falseNegativeWeight = obj.InsidePriorityRatio;
                end

                % C(i,j) is the cost of error: i=true class, j=predicted class
                % Classes: [false,true]
                % matrix: [true negative, false positive; false negative, true positive]
                trainingErrorCost = [0 1;falseNegativeWeight 0];

                [~,svmTrainingOptions] = merge_name_value_pair_argument(...
                    defaultSvmTrainingOptions,...
                    obj.SvmTrainingOptions);

                obj.Svm = fitcsvm(designSample,isInside,...
                    'ClassNames',[false,true],...
                    'Cost',trainingErrorCost,...
                    svmTrainingOptions{:}); 
            end
        end
        
        function [isInside, score] = is_in_candidate_space(obj,designSample)
        %IS_IN_CANDIDATE_SPACE Verification if given design samples are inside
        %   IS_IN_CANDIDATE_SPACE uses the currently defined candidate space to 
        %   determine if given design sample points are inside or outside the candidate 
        %   space.
        %
        %   ISINSIDE = OBJ.IS_IN_CANDIDATE_SPACE(DESIGNSAMPLE) receives the design
        %   samples in DESIGNSAMPLE and returns whether or not they are inside the 
        %   candidate space in ISINSIDE. For ISINSIDE values of 'true', it means the 
        %   respective design is inside the candidate space, while 'false' means it is 
        %   outside.
        %
        %   [ISINSIDE,SCORE] = OBJ.IS_IN_CANDIDATE_SPACE(...) also returns a SCORE value
        %   for each sample point; negative values of SCORE indicate the design sample 
        %   is inside the candidate space, and positive values indicate it is outside. 
        %   Designs with lower/higher SCORE are further from the boundary, with 0 
        %   representing that they are exactly at the boundary.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceSvm
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - ISINSIDE : (nSample,1) logical
        %       - SCORE : (nSample,1) double
        %   
        %   See also predict.

            [isInside,score] = predict(obj.Svm,designSample);
    
            % convert to true (1) / false (0)
            isInside = logical(isInside); 

            % only pick score of how negative designs are (negative score = is positive)
            score = score(:,1);

            % Find Box of Samples Used
            boundingBox = design_bounding_box(obj.DesignSampleDefinition);

            % Assign Value of data outside training set to negative
            outsideTrainingData = ~is_in_design_box(designSample,boundingBox);
            isInside(outsideTrainingData) = false;
            score(outsideTrainingData) = abs(score(outsideTrainingData));
        end

        function obj = grow_candidate_space(obj,growthRate)
        %GROW_CANDIDATE_SPACE Expansion of candidate space by given factor
        %   GROW_CANDIDATE_SPACE will grow the region considered inside the current 
        %   candidate space by the factor given. Said growth is done in a fixed rate 
        %   defined by the input relative to the design space.
        %
        %   OBJ = OBJ.GROW_CANDIDATE_SPACE(GROWTHRATE) will growth the candidate space 
        %   defined in OBJ by a factor of GROWTHRATE. This is an isotropic expansion of 
        %   the candidate space by a factor of the growth rate times the size of the 
        %   design space.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceBase
        %       - GROWTHRATE : double
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceBase
        %   
        %   See also is_in_candidate_space.

            [designSampleExpanded,labelExpanded] = grow_sample_region_positive_label(...
                obj.DesignSampleDefinition,...
                obj.IsInsideDefinition,...
                obj.DesignSpaceLowerBound,...
                obj.DesignSpaceUpperBound,...
                growthRate,...
                obj.GrowthDistanceOptions{:});
            
            % train new candidate space
            obj = obj.define_candidate_space(designSampleExpanded,labelExpanded);
        end

        function plotHandle = plot_candidate_space(obj,figureHandle,varargin)
        %PLOT_CANDIDATE_SPACE Visualization of the boundary of the canidate space 2D/3D
        %   PLOT_CANDIDATE_SPACE allows for the visualization of the boundary of the
        %   candidate space in the given figure. 
        %
        %   PLOTHANDLE = OBJ.PLOT_CANDIDATE_SPACE(FIGUREHANDLE) plots the boundary of
        %   the candidate space in figure FIGUREHANDLE, returning the handle of the 
        %   object plot PLOTHANDLE.
        %
        %   PLOTHANDLE = OBJ.PLOT_CANDIDATE_SPACE(...,NAME,VALUE) allows the 
        %   specification for additional options in the process. For 2D candidate 
        %   spaces, these options should refer to 'contour', and for 3D spaces, they
        %   should refer to 'patch'.
        %
        %   Input:
        %       - OBJ : CandidateSpaceBase
        %       - FIGUREHANDLE : Figure
        %
        %   Output:
        %       - PLOTHANDLE : contour object OR patch object
        %
        %   See also plot_svm_decision_boundary_2d, plot_svm_decision_boundary_3d.

            nDimension = size(obj.ActiveDesign,2);
            if(nDimension==2)
                plotHandle = plot_svm_decision_boundary_2d(figureHandle,obj.Svm,varargin{:});
            elseif(nDimension==3)
                plotHandle = plot_svm_decision_boundary_3d(figureHandle,obj.Svm,varargin{:});
            else
                plotHandle = [];
            end
        end

        function designSample = get.DesignSampleDefinition(obj)
            if(isempty(obj.Svm))
                designSample = [];
            else
                designSample = obj.Svm.X;
            end
        end

        function isInside = get.IsInsideDefinition(obj)
            if(isempty(obj.Svm))
                isInside = [];
            else
                isInside = obj.Svm.Y;
            end
        end
        
        function volume = get.Measure(obj)
            nSample = size(obj.DesignSampleDefinition,1);
            samplingBox = obj.SamplingBox;
            
            volumeSample = sampling_latin_hypercube(samplingBox,nSample);
            isInside = obj.is_in_candidate_space(volumeSample);
            volumeFactor = sum(isInside) / size(isInside,1);
            volume = volumeFactor * prod(samplingBox(2,:) - samplingBox(1,:));
        end

        function isShapeDefinition = get.IsShapeDefinition(obj)
            isShapeDefinition = obj.Svm.IsSupportVector;
        end
    end
end