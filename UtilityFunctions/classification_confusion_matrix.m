function [falsePositive,falseNegative,truePositive,trueNegative] = classification_confusion_matrix(trueLabel,predictedLabel)
%CLASSIFICATION_CONFUSION_MATRIX Confusion Matrix elements for Classification
%   CLASSIFICATION_CONFUSION_MATRIX returns the elements of the classification 
%   confusion matrix for a given classification task with the correct labels  
%   given for benchmarking.
%   These elements are:
%       - False Positives (Type I Error): true label is negative, was classified 
%       as positive
%       - False Negatives (Type II Error): true label is positive, was 
%       classified as negative
%       - True Positives: true label is positive, was classified as positive
%       - True Negatives: true label is negative, was classified as negative
%
%   FALSEPOSITIVE = CLASSIFICATION_CONFUSION_MATRIX(TRUELABEL,PREDICTEDLABEL)
%   receives the correct label of each sample in TRUELABEL and the predicted
%   labels in PREDICTEDLABEL, both of which are given as logical arrays
%   (negative: 'false' / positive: 'true'). The function then returns a logical 
%   flag indicating if a Type I error was commited for each sample in  
%   FALSEPOSITIVE (true if commmited, false if not).
%
%   [FALSEPOSITIVE,FALSENEGATIVE] = CLASSIFICATION_CONFUSION_MATRIX(...) also
%   returns a logical flag indicating if a Type II error was commited for each
%   sample in FALSENEGATIVE (true if commmited, false if not).
%
%   [FALSEPOSITIVE,FALSENEGATIVE,TRUEPOSITIVE] = 
%   CLASSIFICATION_CONFUSION_MATRIX(...) also returns a logical flag indicating 
%   if both correct and predicted labels match as positive for each sample in 
%   TRUEPOSITIVE (true if correct, false if not).
%
%   [FALSEPOSITIVE,FALSENEGATIVE,TRUEPOSITIVE,TRUENEGATIVE] = 
%   CLASSIFICATION_CONFUSION_MATRIX(...) also returns a logical flag indicating 
%   if both correct and predicted labels match as negative for each sample in 
%   TRUENEGATIVE (true if correct, false if not).
%
%   Inputs:
%       - TRUELABEL : (nSample,nTest) logical
%       - PREDICTEDLABEL : (nSample,nTest) logical
%
%   Outputs:
%       - FALSEPOSITIVE : (nSample,nTest) logical
%       - FALSENEGATIVE : (nSample,nTest) logical
%       - TRUEPOSITIVE : (nSample,nTest) logical
%       - TRUENEGATIVE : (nSample,nTest) logical
%
%   See also classification_evaluation_metrics.
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

    % type I error (false positive) - is negative, classified positive
    falsePositive  = (~trueLabel & predictedLabel); 
    % type II error (false negative) - is positive, classified negative
    falseNegative = (trueLabel & ~predictedLabel); 
    % true positive - is positive, classified positive
    truePositive  = (trueLabel & predictedLabel); 
    % true negative - is negative, classified negative
    trueNegative = (~trueLabel & ~predictedLabel); 
end