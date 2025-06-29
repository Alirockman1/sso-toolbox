classdef DesignFastForwardSvm < DesignFastForwardBase
%DESIGNFASTFORWARDSVM Surrogate modeling with Support Vector Machine
%	DESIGNFASTFORWARDSVM uses a Support Vector Machine to create a surrogate
%	model of given labels for a problem. It can predict a score that tells
%	how likely that prediction is to be true, but cannot produce individual
%	deficits for each measure.
%
%   This class is derived from 'DesignFastForwardBase'.
%
%	DESIGNFASTFORWARDSVM specific properties:
%		- SvmModel : 'ClassificationSVM' model used.
%		- SvmTrainingOptions : options to be used when training the SVM model.
%		- SvmBalanceCost : function on how to balance the costs of an given dataset.
%
%	See also DesignFastForwardBase, DesignFastForwardAnn.
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

	properties
		%SVMMODEL Support Vector Machine
        %   SVMMODEL is the object that contains the support vector machine used for
        %   the prediction of designs.
        %
        %   SVMMODEL : ClassificationSVM object
        %
        %   See also ClassificationSVM.
		SvmModel

		%SVMTRAININGOPTIONS Options for SVM training 
        %   SVMTRAININGOPTIONS is a cell that contains further options for the
        %   training of the support vector machine. In its default value, it
        %   contains the following parameters:
        %       - Standardize : true
        %       - KernelFunction : 'gaussian'
        %       - Verbose : 0
        %   These and other options of 'fitcsvm' can be chosen as desired, with
        %   the following two exceptions:
        %       - 'Cost' should not be changed, as it is already defined (see 
        %       'SvmBalanceCost' property).
        %       - 'ClassNames', which is fixed to [false,true].
        %   SVMTRAININGOPTIONS should be defined as a cell array with name-value pairs 
        %   or a structure.
        %
        %   SVMTRAININGOPTIONS : cell with name-value pairs OR structure
        %
        %   See also fitcsvm, SvmBalanceCost, ClassificationSVM.
		SvmTrainingOptions

		%SVMBALANCECOST Balance false positive/negative cost in SVM training
        %   SVMBALANCECOST determines the weight of making false positive errors
        %   when training the Support Vector Machine that defines the candidate
        %   space. The weight of making false negative errors is fixed at 1.
        %   SVMBALANCECOST can be determined as either a fixed number or as a 
        %   function which takes the number of positive/negatives points in the
        %   training dataset. In the second case, the weight should increase
        %   as the number of positive samples becomes higher than the number
        %   of negative samples.
        %
        %   SVMBALANCECOST : double OR function handle with the format
        %   'weight = f(nPositiveTraining,nNegativeTraining)' 
        %
        %   See also fitcsvm, SVMTRAININGOPTIONS, ClassificationSVM.
		SvmBalanceCost
	end

	properties(SetAccess = protected)
		%MEASURETRAIN Measures of design samples used in the training of current model
        %   MEASURETRAIN are the measures of the design samples used in the training of 
        %	the current model. This can be either performance measures or physical 
        %	feasibility measures.
        %
        %   MEASURETRAIN : (nSampleTrain,nMeasure) double
        %
        %   See also DesignSampleTrain, LabelTrain.
		MeasureTrain
	end

	properties(SetAccess = protected, Dependent)
		%DESIGNSAMPLETRAIN Design samples used in the training of the current model
        %   DESIGNSAMPLETRAIN are the design samples used in the training of the current
        %	model.
        %
        %   DESIGNSAMPLETRAIN : (nSampleTrain,nDesignVariable) double
        %
        %   See also MeasureTrain, LabelTrain.
		DesignSampleTrain
	end

	methods
		function obj = DesignFastForwardSvm(measureLowerLimit,measureUpperLimit,varargin)
		%DESIGNFASTFORWARDSVM Constructor
		%	DESIGNFASTFORWARDSVM initializes a design fast-forward model with a support
		%	vector machine being used for predictions.
		%
		%	OBJ = DESIGNFASTFORWARDSVM(MEASURELOWERLIMIT,MEASUREUPPERLIMIT) creates an
		%	initial object with the critical lower/upper limits for the measures 
		%	MEASURELOWERLIMIT and MEASUREUPPERLIMIT. Options for creating the net and
		%	for training are set to default. The model and training information are 
		%	initialized as empty.
		%
		%	OBJ = DESIGNFASTFORWARDSVM(MEASURELOWERLIMIT,MEASUREUPPERLIMIT,NAME1,
		%	VALUE1,...) allows one to also set additional options for the creation of 
		%	the net and its training.
		%		- SvmTrainingOptions : options for the training of the SVM; see more
        %   	on the respective property header.
        %   	- SvmBalanceCost : weight factor definition for false positive 
        %   	errors comitted during training. See more on the respective property
        %   	header.
		%
		%	Inputs:
		%		- MEASURELOWERLIMIT : (1,nMeasure) double
		%		- MEASUREUPPERLIMIT : (1,nMeasure) double
		%
		%	Output:
		%		- OBJ : DesignFastForwardSvm object
		%
		%	See also ClassificationSVM, SvmTrainingOptions, SvmBalanceCost.
			obj.MeasureLowerLimit = measureLowerLimit;
			obj.MeasureUpperLimit = measureUpperLimit;
			obj.SvmModel = [];

			parser=inputParser;
			parser.addParameter('SvmTrainingOptions',{},@(x)iscell(x));
			parser.addParameter('SvmBalanceCost',...
				@(nPositive,nNegative)[min([max([nPositive/nNegative,0.1]),10])],...
                @(x)isa(x,'function_handle'));
			parser.parse(varargin{:});

			defaultSvmTrainingOptions = {...
				'Standardize',true,...
				'KernelFunction','gaussian',...
				'OptimizeHyperparameters','auto',...
				'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus',...
			        'UseParallel',true,'MaxObjectiveEvaluations',50,'Verbose',0),...
				'Verbose',0,...
				'NumPrint',10};
			[~,obj.SvmTrainingOptions] = merge_name_value_pair_argument(defaultSvmTrainingOptions,parser.Results.SvmTrainingOptions);

			obj.SvmBalanceCost = parser.Results.SvmBalanceCost;
		end

		function obj = train_model(obj,designSample,measure)
		%TRAIN_MODEL Train design fast-forward model anew with given samples
		%	TRAIN_MODEL initiates the model object by training it with the new samples
		%	given. If there was a model before, it gets overwritten.
		%
		%	OBJ = OBJ.TRAIN_MODEL(DESIGNSAMPLE,MEASURE) uses the design samples 
		%	DESIGNSAMPLE and their respective measures MEASURE to train a new model
		%	and return the new object OBJ.
		%
		%	Inputs:
		%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
		%		- MEASURE : (nSample,nMeasure) double
		%
		%	Output:
		%		- OBJ : DesignFastForward object
		%
		%	See also update_model.

			% associate two classes to data: +1 for bad, -1 for good
			obj.MeasureTrain = measure;
			label = design_fulfills_limit_criteria(measure,obj.MeasureLowerLimit,obj.MeasureUpperLimit);

			% Weighting factor of Type I error (false positive - is bad design, predicted good design)
			nNegativeSample = sum(~label);
			nPositiveSample = sum(label);
			if(isa(obj.SvmBalanceCost,'function_handle'))
                falseNegativeWeight = obj.SvmBalanceCost(nPositiveSample,nNegativeSample);
            else
                falseNegativeWeight = obj.SvmBalanceCost;
            end

            % C(i,j) is the cost of error: i=true class, j=predicted class
            trainingErrorCost = [0 1;falseNegativeWeight 0];

			%% SVM Training
			if(all(label == label(1))) % single class
			    obj.SvmModel = fitcsvm(designSample,label,'Standardize',true,'KernelFunction','gaussian'); 
			else % two classes
			    obj.SvmModel = fitcsvm(designSample,label,'Cost',trainingErrorCost,'ClassNames',[false,true],...
			       	obj.SvmTrainingOptions{:});
			    gcf; close; gcf; close;
			end
		end

		function obj = update_model(obj,newDesignSample,newMeasure)
		%UPDATE_MODEL Update design fast-forward model with new given samples
		%	UPDATE_MODEL updates the model object by training it with the new samples
		%	given. The previously used samples to train the current model are equally
		%	considered.
		%
		%	OBJ = OBJ.UPDATE_MODEL(NEWDESIGNSAMPLE,NEWMEASURE) uses the design samples 
		%	NEWDESIGNSAMPLE and their respective measures NEWMEASURE to update the  
		%	current model and return the new object OBJ.
		%
		%	Inputs:
		%		- NEWDESIGNSAMPLE : (nSample,nDesignVariable) double
		%		- NEWMEASURE : (nSample,nMeasure) double
		%
		%	Output:
		%		- OBJ : DesignFastForward object
		%
		%	See also train_model.

			trainDesignSample = [obj.DesignSampleTrain;newDesignSample];
			trainMeasure = [obj.MeasureTrain;newMeasure];
			obj = obj.train_model(trainDesignSample,trainMeasure);
		end

		function [score,output] = predict_designs(obj,designSample)
		%PREDICT_DESIGNS Predict deficit of given design samples with the current model
		%	PREDICT_DESIGNS uses the current model to predict the measure deficit of the 
		%	given design samples. 
		%
		%	DEFICIT = OBJ.PREDICT_DESIGNS(DESIGNSAMPLE) predicts for the design samples
		%	DESIGNSAMPLE the associated measure deficit DEFICIT.
		%
		%	[DEFICIT,OUTPUT] = OBJ.PREDICT_DESIGNS(DESIGNSAMPLE) also returns the output
		%	of said prediction; here, it returns the predicted SVM labels.
		%
		%	Inputs:
		%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
		%
		%	Output:
		%		- DEFICIT : (nSample,nMeasure) double
		%		- OUTPUT : class-specific
		%
		%	See also .

			[output.PredictedLabel,score] = predict(obj.SvmModel,designSample);
			
			% only pick score of how "bad" designs are (negative values = is good)
			score = score(:,1);
		end

		function DesignSampleTrain = get.DesignSampleTrain(obj)
			DesignSampleTrain = obj.SvmModel.X;
		end
	end
end