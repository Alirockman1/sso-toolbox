classdef CandidateSpaceDecisionTree < CandidateSpaceBase
%CANDIDATESPACEDECISIONTREE Candidate Space defined by a Support Vector Machine
%   CANDIDATESPACEDECISIONTREE can be used to define candidate spaces for  
%   component solution space computation. In this case, the candidate space is 
%   defined by a decision tree trained on the labeled data given for design 
%   sample points.
%
%   This class is derived from CandidateSpaceBase.
%
%   CANDIDATESPACEDECISIONTREE properties:
%       - DesignSampleDefinition : design sample points used in the candidate 
%       space definition.
%       - IsInsideDefinition : logical labels of the design sample points used   
%       in the candidate space definition regarding whether they are inside or 
%       outside the candidate space (true = inside, false = outside).
%       - IsShapeDefinition : logical labels of the design sample points used in 
%       the candidate space definition regarding whether these designs are 
%       directly related to the shape of the candidate space (or in other words,
%       the design sample points that define the boundary itself and are used
%       to identify if a design is inside/outside this space).
%       - ActiveDesign : sample points that are labeled as inside.
%       - Measure : measure of the candidate space (usually area/volume/etc).
%       - SamplingBox : bounding box around the internal region of the candidate
%       space that can be used to better attempt to sample inside said region.
%       - SamplingBoxSlack : a value used to give the defined amount of slack
%       for the SamplingBox, varying between the strict bounding box around the 
%       internal region and a larger bounding box which may contain negative
%       designs in its edges. 
%       - GrowthDistanceOptions : options to be used when computing the distance 
%       that the candidate space will grow to. In its default base 
%       implementation, these are options to 'knnsearch'.
%       - DecisionTree : 
%       - DecisionTreeTrainingOptions : 
%       - Boundary
%
%   CANDIDATESPACEDECISIONTREE methods:
%       - define_candidate_space : create a candidate space based on design 
%       samples that are labeled as inside/outside.
%       - grow_candidate_space : expand the candidate space by a given factor.
%       - is_in_candidate_space : verify if given design samples are inside 
%       the candidate space.
%
%   See also ClassificationTree.
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
        %   For CandidateSpaceDecisionTree, MEASURE is the volume of the inside region 
        %   which is estimated with the Monte Carlo integration method.
        %
        %   MEASURE : double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition. 
        Measure
    end

    properties (SetAccess = protected)
        %DECISIONTREE Classification Decision Tree
        %   DECISIONTREE is the object that contains the decision tree classifier used 
        %   for the definition and classification of the candidate space.
        %
        %   DECISIONTREE : ClassificationTree 
        %
        %   See also fitctree, ClassificationTree.
        DecisionTree

        %DECISIONTREETRAININGOPTIONS Options for Tree training during definition
        %   DECISIONTREETRAININGOPTIONS is a cell that contains further options for the
        %   training of the Decision Tree. In its default value, it
        %   contains the following parameters:
        %       - 'ScoreTransform' : 'symmetric'
        %   Any options from 'fitctree' can be chosen as desired, with the following two
        %   exceptions:
        %       - 'Cost' should not be changed, as it is already defined (see 
        %       'InsidePriorityRatio' property).
        %       - 'ClassNames', which is fixed to [false,true].
        %   DECISIONTREETRAININGOPTIONS should be defined as a cell with name-value pair 
        %   arguments or a structure.
        %
        %   DECISIONTREETRAININGOPTIONS : cell with name-value pairs OR structure
        %
        %   See also fitctree, InsidePriorityRatio.
        DecisionTreeTrainingOptions

        %INSIDEPRIORITYRATIO Cost of false negative classification in Tree training
        %   INSIDEPRIORITYRATIO determines the weight of making false negative errors
        %   (is inside, classified outside) when training the Decision Tree that defines
        %   the candidate space. The weight of making false positive errors (is outside, 
        %   classified inside) is fixed at 1.
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
        %   See also fitctree, DecisionTreeTrainingOptions.
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
        function obj = CandidateSpaceDecisionTree(designSpaceLowerBound,designSpaceUpperBound,varargin)
        %CandidateSpaceDecisionTree Constructor
        %   CandidateSpaceDecisionTree is a constructor initializes an object of this 
        %   class.
        %
        %   OBJ = CandidateSpaceDecisionTree(DESIGNSPACELOWERBOUND,
        %   DESIGNSPACEUPPERBOUND) creates an object of the CandidateSpaceDecisionTree
        %   class and sets its design space boundaries to its input values 
        %   DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND. Other properties are set to
        %   empty.
        %
        %   OBJ = CandidateSpaceDecisionTree(...,NAME,VALUE,...) also allows one to set
        %   specific options for the object. This can be 
        %       - 'SamplingBoxSlack' : where the boundaries of the sampling box will be
        %       relative to the strictest bounding box and the most relaxed bounding 
        %       box. A value of 0 means no slack and therefore the sampling box will be
        %       the most strict one possible, and 1 means the sampling box will be the 
        %       most relaxed.
        %       - 'DecisionTreeTrainingOptions' : options for the training of the 
        %       decision tree; see more on the respective property header.
        %       - 'InsidePriorityRatio' : weight factor definition for false negative 
        %       errors comitted during training. See more on the respective property
        %       header.
        %
        %   Inputs:
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %       - 'SamplingBoxSlack' : double
        %       - 'DecisionTreeTrainingOptions' : cell with name-value pairs OR 
        %       structure
        %       - 'InsidePriorityRatio' : function handle OR double
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceDecisionTree
        %   
        %   See also DecisionTreeTrainingOptions, InsidePriorityRatio.

            % parse inputs
            parser = inputParser;
            parser.addRequired('designSpaceLowerBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addRequired('designSpaceUpperBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addParameter('SamplingBoxSlack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
            parser.addParameter('DecisionTreeTrainingOptions',{},@(x)iscell(x));
            parser.addParameter('InsidePriorityRatio',1,@(x)isa(x,'function_handle')||isnumeric(x));
            parser.addParameter('GrowthDistanceOptions',{},@(x)iscell(x));
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            % initialize properties
            obj.DecisionTree = [];
            obj.DesignSpaceLowerBound = parser.Results.designSpaceLowerBound;
            obj.DesignSpaceUpperBound = parser.Results.designSpaceUpperBound;
            obj.InsidePriorityRatio = parser.Results.InsidePriorityRatio;
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;
            obj.DecisionTreeTrainingOptions = parser.Results.DecisionTreeTrainingOptions;

            defaultGrowthDistanceOptions = {};
            [~,obj.GrowthDistanceOptions] = merge_name_value_pair_argument(...
                defaultGrowthDistanceOptions,...
                parser.Results.GrowthDistanceOptions);
        end
        
        function obj = define_candidate_space(obj,designSample,isInside)
        %DEFINE_CANDIDATE_SPACE Initial definition of the candidate space
        %   DEFINE_CANDIDATE_SPACE uses labeled design samples to define the inside / 
        %   outside regions of the candidate space.
        %   For CandidateSpaceDecisionTree, this means the labeled design samples are 
        %   used to train a Decision Tree.
        %
        %   OBJ = OBJ.DEFINE_CANDIDATE_SPACE(DESIGNSAMPLE) receives the design samle
        %   points in DESIGNSAMPLE and returns a candidate space object OBJ with the new
        %   definition, assuming all designs are inside the candidate space.
        %
        %   OBJ = OBJ.DEFINE_CANDIDATE_SPACE(DESIGNSAMPLE,ISINSIDE) additionally 
        %   receives the inside/outside (true/false) labels of each design point in 
        %   ISINSIDE.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceDecisionTree
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceDecisionTree
        %   
        %   See also fitctree, is_in_candidate_space.

            if(nargin<3)
                isInside = true(size(designSample,1),1);
            end

            defaultDecisionTreeTrainingOptions = {...
                ...'HyperparameterOptimizationOptions',struct('Verbose',0),...
                ...'OptimizeHyperparameters','auto',...
                'ScoreTransform','symmetric'};

            % Decision Tree Training
            if(all(isInside == isInside(1))) % single class
                obj.DecisionTree = fitctree(designSample,isInside,defaultDecisionTreeTrainingOptions{:}); 
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

                [~,decisionTreeTrainingOptions] = merge_name_value_pair_argument(...
                    defaultDecisionTreeTrainingOptions,...
                    obj.DecisionTreeTrainingOptions);

                obj.DecisionTree = fitctree(designSample,isInside,...
                    'ClassNames',[false,true],...
                    'Cost',trainingErrorCost,...
                    decisionTreeTrainingOptions{:}); 
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
        %       - OBJ : CandidateSpaceDecisionTree
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - ISINSIDE : (nSample,1) logical
        %       - SCORE : (nSample,1) double
        %   
        %   See also predict.

            [isInside,score] = predict(obj.DecisionTree,designSample);
    
            % convert to true (1) / false (0)
            isInside = logical(isInside); 

            % negative score = is inside
            score = score(:,1);
            score(isInside) = -abs(score(isInside));
            score(~isInside) = abs(score(~isInside));

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

        function designSample = get.DesignSampleDefinition(obj)
            if(isempty(obj.DecisionTree))
                designSample = [];
            else
                designSample = obj.DecisionTree.X;
            end
        end

        function isInside = get.IsInsideDefinition(obj)
            if(isempty(obj.DecisionTree))
                isInside = [];
            else
                isInside = obj.DecisionTree.Y;
            end
        end

        function isShapeDefinition = get.IsShapeDefinition(obj)
            if(all(obj.IsInsideDefinition(1)==obj.IsInsideDefinition))
                isShapeDefinition = true(size(obj.IsInsideDefinition,1),1);
            else
                isInBoundary = design_find_boundary_samples(obj.DesignSampleDefinition,obj.IsInsideDefinition);

                [~,iMinActiveDesignVariable] = min(obj.ActiveDesign,[],1);
                [~,iMaxActiveDesignVariable] = max(obj.ActiveDesign,[],1);
                nSample = size(obj.DesignSampleDefinition,1);
                globalIndex = convert_index_base(...
                    obj.IsInsideDefinition,...
                    [iMinActiveDesignVariable';iMaxActiveDesignVariable'],...
                    'backward');
                isMaxMin = ismember((1:nSample)',globalIndex);

                isShapeDefinition = isInBoundary | isMaxMin;
            end
        end
        
        function volume = get.Measure(obj)
            nSample = size(obj.DesignSampleDefinition,1);
            samplingBox = obj.SamplingBox;
            
            % use Monte Carlo approximation
            volumeSample = sampling_latin_hypercube(samplingBox,nSample);
            isInside = obj.is_in_candidate_space(volumeSample);
            volumeFactor = sum(isInside) / size(isInside,1);
            volume = volumeFactor * prod(samplingBox(2,:) - samplingBox(1,:));
        end
    end
end