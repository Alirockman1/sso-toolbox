classdef (Abstract) DesignFastForwardBase
%DESIGNFASTFORWARDBASE Design Fast-Forward abstract class for surrogate modeling
%	DESIGNFASTFORWARDBASE provides the interfaces and basic functions necessary 
%	for an implementation of the DesignFastForward-type class. These classes and
%	objects can then be used for surrogate modeling of designs in the context
%	of top-down design, predicting whether the given designs are good/physically
%	feasible or bad/physically infeasible.
%
%   As an abstract class, particular implementations have to be created using 
%   this as base for one to be able to create objects.
%
%	DESIGNFASTFORWARDBASE properties:
%		- MeasureLowerLimit : critical lower limit for the considered measures.
%		- MeasureUpperLimit : critical upper limit for the considered measures.
%		- DesignSampleTrain : design samples used for training the current 
%		model.
%		- MeasureTrain : measures used for training the current model.
%		- LabelTrain : label for the design samples used for training the
%		current model.
%
%	DESIGNFASTFORWARDBASE methods:
%		- train_model : train a new model, discarting any model built before.
%		- update_model : update an existing model with new data.
%		- predict_designs : predict the deficit of the given design samples.
%		- set_critical_limit : set the critical lower/upper limits for the 
%		measures.
%		- classify_measure : given measures, find their label/score/deficit.
%
%	See also DesignFastForwardAnn, DesignFastForwardSvm.
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
		%MEASURELOWERLIMIT Critical lower limit for the measures
        %   MEASURELOWERLIMIT is the critical lower limit for the measures.
        %
        %   MEASURELOWERLIMIT : (1,nMeasure) double
        %
        %   See also MeasureUpperLimit.
		MeasureLowerLimit

		%MEASUREUPPERLIMIT Critical upper limit for the measures
        %   MEASUREUPPERLIMIT is the critical upper limit for the measures.
        %
        %   MEASUREUPPERLIMIT : (1,nMeasure) double
        %
        %   See also MeasureLowerLimit.
		MeasureUpperLimit
	end

	properties (SetAccess = protected, Abstract)
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
	end

	properties (SetAccess = protected, Dependent)
		%LABELTRAIN Labels of design samples used in the training of current model
        %   LABELTRAIN are the labels of the design samples used in the training of the
        %	current model. Its value is 'true' for positive samples (those that have 
        %	measures within the specified limits) and 'false' for negative samples 
        %	(those that have measures outside the specified limits).
        %
        %   LABELTRAIN : (nSampleTrain,1) logical
        %
        %   See also DesignSampleTrain, MeasureTrain.
		LabelTrain

		%
		BoundingBoxPositive
	end

	methods (Abstract)
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
		obj = train_model(obj,designSample,measure)
		
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
		obj = update_model(obj,newDesignSample,newMeasure)
		
		%PREDICT_DESIGNS Predict deficit of given design samples with the current model
		%	PREDICT_DESIGNS uses the current model to predict the measure deficit of the 
		%	given design samples. 
		%
		%	DEFICIT = OBJ.PREDICT_DESIGNS(DESIGNSAMPLE) predicts for the design samples
		%	DESIGNSAMPLE the associated measure deficit DEFICIT.
		%
		%	[DEFICIT,OUTPUT] = OBJ.PREDICT_DESIGNS(DESIGNSAMPLE) also returns the output
		%	of said prediction; how this is handled is defined on a class-by-class 
		%	basis.
		%
		%	Inputs:
		%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
		%
		%	Output:
		%		- DEFICIT : (nSample,nMeasure) double
		%		- OUTPUT : class-specific
		%
		%	See also .
		[deficit,output] = predict_designs(obj,designSample)
	end

	methods
		function obj = set_critical_limit(obj,measureLowerLimit,measureUpperLimit)
		%SET_CRITICAL_LIMIT Set the lower/upper critical limits for the measures
		%	SET_CRITICAL_LIMIT allows one to set up the lower/upper critical limit 
		%	values for the considered measures.
		%
		%	OBJ = OBJ.SET_CRITICAL_LIMIT(MEASURELOWERLIMIT,MEASUREUPPERLIMIT)
		%	sets the lower/upper limits of the  measures (MEASURELOWERLIMIT and 
		%	MEASUREUPPERLIMIT) for the DesignFastForward object OBJ. 
		%
		%	Inputs:
		%		- MEASURELOWERLIMIT : (1,nMeasure) double
		%		- MEASUREUPPERLIMIT : (1,nMeasure) double
		%
		%	Output:
		%		- OBJ : DesignFastForward
		%
		%	See also .
			if(size(measureLowerLimit,1)~=1 || size(measureUpperLimit,1)~=1 || size(measureLowerLimit,2)~=size(measureUpperLimit,2))
				error('DesignFastForward:CriticalLimitSizeError','Incompatible sizes for critical limits of the measures');
			end

			obj.MeasureLowerLimit = measureLowerLimit;
			obj.MeasureUpperLimit = measureUpperLimit;
		end

		function [label,score,deficit] = classify_measure(obj,measure)
		%CLASSIFY_MEASURE Classify given measures with the model's current limits
		%	CLASSIFY_MEASURE uses the current critical lower/upper limits of the model
		%	to classify the given measures.
		%
		%	LABEL = OBJ.CLASSIFY_MEASURE(MEASURE) gives the positive/negative ('true'/
		%	'false') label LABEL for each design sample measure MEASURE.
		%
		%	[LABEL,SCORE] = OBJ.CLASSIFY_MEASURE(...) also returns the SCORE for each
		%	design sample.
		%
		%	[LABEL,SCORE,DEFICIT] = OBJ.CLASSIFY_MEASURE(...) also returns the DEFICIT
		%	for each measure.
		%
		%	Input:
		%		- MEASURE : (nSample,nMeasure) double
		%
		%	Output:
		%		- LABEL : (nSample,1) logical
		%		- SCORE : (nSample,1) double
		%		- DEFICIT : (nSample,nMeasure) double
		%
		%	See also design_measure_to_deficit, design_deficit_to_label_score.
			deficit = design_measure_to_deficit(measure,obj.MeasureLowerLimit,obj.MeasureUpperLimit);
			[label,score] = design_deficit_to_label_score(deficit);
		end

		function boundingBox = get.BoundingBoxPositive(obj)
			boundingBox = design_bounding_box(obj.DesignSampleTrain,obj.LabelTrain);
		end

		function labelTrain = get.LabelTrain(obj)
			labelTrain = obj.classify_measure(obj.MeasureTrain);
		end
	end
end