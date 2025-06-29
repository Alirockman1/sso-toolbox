function [accuracy,precision,recall,specificity,f1score] = classification_evaluation_metrics(falsePositive,falseNegative,truePositive,trueNegative)
%CLASSIFICATION_EVALUATION_METRICS Important metrics of classification problem
%   CLASSIFICATION_EVALUATION_METRICS returns a handful of common evaluation 
%   metrics for a given classification task. Given the number (or individual 
%   arrays) of false positives / false negatives / true positives / true 
%   negatives in the task, CLASSIFICATION_EVALUATION_METRICS can compute the 
%   most famous metrics used to evaluate the quality of the classifier.
%
%   ACCURACY = CLASSIFICATION_EVALUATION_METRICS(FALSEPOSITIVE,FALSENEGATIVE,
%   TRUEPOSITIVE,TRUENEGATIVE) receives the elements of the confusion matrix
%   as either a given total number or as boolean arrays (where 'true' means that
%   sample point is of that element) in FALSEPOSITIVE, FALSENEGATIVE, 
%   TRUEPOSITIVE, and TRUENEGATIVE. It returns the ACCURACY of the 
%   classification, which is the ratio of the number of true classifications by 
%   the total number of samples (true / all).
%
%   [ACCURACY,PRECISION] = CLASSIFICATION_EVALUATION_METRICS(...) also returns
%   the PRECISION of the classification, which is the ratio of the number of 
%   true positive classifications by the number of samples that got the positive
%   label (true positive / predicted positive).
%
%   [ACCURACY,PRECISION,RECALL] = CLASSIFICATION_EVALUATION_METRICS(...) also 
%   returns the RECALL of the classification, which is the ratio of the number  
%   of true positive classifications by the number of samples that are the 
%   positive label (true positive / is positive).
%
%   [ACCURACY,PRECISION,RECALL,SPECIFICITY] = 
%   CLASSIFICATION_EVALUATION_METRICS(...) also returns the SPECIFICITY of the 
%   classification, which is the ratio of the number of true negative 
%   classification by the number of samples that are the negative label
%   (true negative / is negative).
%
%   [ACCURACY,PRECISION,RECALL,SPECIFICITY,F1SCORE] = 
%   CLASSIFICATION_EVALUATION_METRICS(...) also returns the F1-score of the 
%   classification in F1SCORE, which is the harmonic mean of precision and 
%   recall.
%
%   Inputs:
%       - FALSEPOSITIVE : (1,nTest) integer OR (nSample,nTest) logical
%       - FALSENEGATIVE : (1,nTest) integer OR (nSample,nTest) logical
%       - TRUEPOSITIVE : (1,nTest) integer OR (nSample,nTest) logical
%       - TRUENEGATIVE : (1,nTest) integer OR (nSample,nTest) logical
%
%   Outputs:
%       - ACCURACY : (1,nTest) double
%       - PRECISION : (1,nTest) double
%       - RECALL : (1,nTest) double
%       - SPECIFICITY : (1,nTest) double
%       - F1SCORE : (1,nTest) double
%
%   See also classification_confusion_matrix.
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

    % get the totals
    nFalsePositive = sum(falsePositive,1);
    nFalseNegative = sum(falseNegative,1);
    nTruePositive = sum(truePositive,1);
    nTrueNegative = sum(trueNegative,1);
    
    % compute evaluation metrics
    accuracy = (nTrueNegative+nTruePositive)./(nFalsePositive+nFalseNegative+nTruePositive+nTrueNegative); % true / all
    precision = nTruePositive./(nTruePositive+nFalsePositive); % true positive / predicted positive
    recall = nTruePositive./(nTruePositive+nFalseNegative); % true positive / is positive
    specificity = nTrueNegative./(nTrueNegative+nFalsePositive); % true negative / is negative
    f1score = harmmean([precision;recall],1);
end