classdef DesignFastForwardAnn < DesignFastForwardBase
%DESIGNFASTFORWARDANN Surrogate modeling with Artificial Neural Networks
%	DESIGNFASTFORWARDANN uses an Aritifical Neural Network to create a surrogate
%	model of given design measures for a problem.
%
%   This class is derived from 'DesignFastForwardBase'.
%
%	DESIGNFASTFORWARDANN specific properties:
%		- AnnModel : 'cascadeforwardnet' model used.
%		- AnnNetOptions : options to be used when creating the net.
%		- AnnTrainingOptions : options to be used when training the model.
%
%	See also DesignFastForwardBase, DesignFastForwardSvm.
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

	properties (SetAccess = protected)
		%ANNMODEL Artificial Neural Network Model
		%	ANNMODEL is a 'cascadeforwardnet' object trained on the given design sample
		%	data.
		%
		%	ANNMODEL : network
		%
		%	See also cascadeforwardnet.
		AnnModel
		
		%DESIGNSAMPLETRAIN Design samples used in the training of the current model
        %   DESIGNSAMPLETRAIN are the design samples used in the training of the current
        %	model.
        %
        %   DESIGNSAMPLETRAIN : (nSampleTrain,nDesignVariable) double
        %
        %   See also MeasureTrain, LabelTrain.
		DesignSampleTrain
		
		%MEASURETRAIN Measures of design samples used in the training of current model
        %   MEASURETRAIN are the measures of the design samples used in the training of 
        %	the current model. This can be either performance measures or physical 
        %	feasibility measures.
        %
        %   MEASURETRAIN : (nSampleTrain,nMeasure) double
        %
        %   See also DesignSampleTrain, LabelTrain.
		MeasureTrain
		
		%ANNNETOPTIONS Options for the net object
		%	ANNNETOPTIONS are the options for the 'cascadeforwardnet' object.
		%
		%	ANNNETOPTIONS : name-value pair cell
		%
		%	See also cascadeforwardnet, AnnTrainingOptions.
		AnnNetOptions
		
		%ANNTRAININGOPTIONS Options for the training operation with the sample data
		%	ANNTRAININGOPTIONS are the options for the 'train' operation with the
		%	Artificial Neural Network model.
		%
		%	ANNTRAININGOPTIONS : name-value pair cell
		%
		% 	See also train.
		AnnTrainingOptions
	end

	methods
		function obj = DesignFastForwardAnn(measureLowerLimit,measureUpperLimit,varargin)
		%DESIGNFASTFORWARDANN Constructor
		%	DESIGNFASTFORWARDANN initializes a design fast-forward model with an 
		%	artificial neural network being used for predictions.
		%
		%	OBJ = DESIGNFASTFORWARDANN(MEASURELOWERLIMIT,MEASUREUPPERLIMIT) creates an
		%	initial object with the critical lower/upper limits for the measures 
		%	MEASURELOWERLIMIT and MEASUREUPPERLIMIT. Options for creating the net and
		%	for training are set to default. The model and training information are 
		%	initialized as empty.
		%
		%	OBJ = DESIGNFASTFORWARDANN(MEASURELOWERLIMIT,MEASUREUPPERLIMIT,NAME1,
		%	VALUE1,...) allows one to also set additional options for the creation of 
		%	the net and its training.
		%		- AnnNetOptions : options for the net object as 'cascadeforwardnet'.
		%		- AnnTrainingOptions : options for the training of the model as 'train'.
		%
		%	Inputs:
		%		- MEASURELOWERLIMIT : (1,nMeasure) double
		%		- MEASUREUPPERLIMIT : (1,nMeasure) double
		%
		%	Output:
		%		- OBJ : DesignFastForwardAnn object
		%
		%	See also cascadeforwardnet, train.
			obj.MeasureLowerLimit = measureLowerLimit;
			obj.MeasureUpperLimit = measureUpperLimit;
			obj.AnnModel = [];
			obj.DesignSampleTrain = [];
			obj.MeasureTrain = [];

			parser=inputParser;
			parser.addParameter('AnnNetOptions',{},@(x)iscell(x));
			parser.addParameter('AnnTrainingOptions',{},@(x)iscell(x));
			parser.parse(varargin{:});

			defaultAnnNetOptions = {};
			[~,obj.AnnNetOptions] = merge_name_value_pair_argument(defaultAnnNetOptions,parser.Results.AnnNetOptions);

			defaultAnnTrainingOptions = {'useParallel','yes','showResources','no'};
			[~,obj.AnnTrainingOptions] = merge_name_value_pair_argument(defaultAnnTrainingOptions,parser.Results.AnnTrainingOptions);
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

			% Find indices where prediction of QoIs is applicable
			isValidXTrain = ~any(isnan(designSample)|isinf(designSample),2);
			isValidYTrain = ~any(isnan(measure)|isinf(measure),2);
			isValidTrain = isValidXTrain & isValidYTrain;
			
			% Create Artificial Neural Network for Quantities of Interest
			net = cascadeforwardnet(obj.AnnNetOptions{:});
			obj.DesignSampleTrain = designSample(isValidTrain,:);
			obj.MeasureTrain = measure(isValidTrain,:);
			obj.AnnModel = train(net,obj.DesignSampleTrain',obj.MeasureTrain',...
                obj.AnnTrainingOptions{:});
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

		function [deficit,output] = predict_designs(obj,designSample)
		%PREDICT_DESIGNS Predict deficit of given design samples with the current model
		%	PREDICT_DESIGNS uses the current model to predict the measure deficit of the 
		%	given design samples. 
		%
		%	DEFICIT = OBJ.PREDICT_DESIGNS(DESIGNSAMPLE) predicts for the design samples
		%	DESIGNSAMPLE the associated measure deficit DEFICIT.
		%
		%	[DEFICIT,OUTPUT] = OBJ.PREDICT_DESIGNS(DESIGNSAMPLE) also returns the output
		%	of said prediction; here, it returns a structure with the estimated 
		%	measures.
		%
		%	Inputs:
		%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
		%
		%	Output:
		%		- DEFICIT : (nSample,nMeasure) double
		%		- OUTPUT : class-specific
		%
		%	See also .

			output.Measure = obj.AnnModel(designSample')';
			deficit = design_measure_to_deficit(output.Measure,obj.MeasureLowerLimit,obj.MeasureUpperLimit);
		end
	end
end