%% Cleanup
fclose all;
close all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(4);


%% Documentation / Archive
RNGstate = rng;
saveFolder = save_diary_files(mfilename);
goldenratio = (1+sqrt(5))/2;
figureSize = [goldenratio 1]*8.5;


%% create separable region
designSpaceLowerBound = [0 0];
designSpaceUpperBound = [1 1];

designSample = designSpaceLowerBound + rand(1000,2).*(designSpaceUpperBound-designSpaceLowerBound);

systemFunction = @distance_to_center;
systemParameter = [0.5 0.5];
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

performanceLowerLimit = 0.15;
performanceUpperLimit = 0.3;
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);
labelSample = design_deficit_to_label_score(designEvaluator.evaluate(designSample));


%% plot
figure;
plot(designSample(labelSample,1),designSample(labelSample,2),'g.');
hold all;
plot(designSample(~labelSample,1),designSample(~labelSample,2),'r.');
grid minor;


%% train candidate space
candidateSpace = CandidateSpaceDelaunay(designSpaceLowerBound,designSpaceUpperBound);
candidateSpace = candidateSpace.define_candidate_space(designSample,labelSample);
isShapeDefinition = candidateSpace.IsShapeDefinition;

figure;
plot(designSample(labelSample,1),designSample(labelSample,2),'g.');
hold all;
plot(designSample(~labelSample,1),designSample(~labelSample,2),'r.');
candidateSpace.plot_candidate_space(gcf,'FaceColor','g','FaceAlpha',0.5,'EdgeColor','k');
plot(designSample(isShapeDefinition,1),designSample(isShapeDefinition,2),'bo');
grid minor;
legend({'Inside Points','Outside Points','Candidate Space Inside Region','Shape Points'});


%% update candidate space
% remove designs close to the border
performanceLowerLimitNew = 0.15;
performanceUpperLimitNew = 0.25;
designEvaluatorNew = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimitNew,performanceUpperLimitNew);
designSampleNew = designSpaceLowerBound + rand(1000,2).*(designSpaceUpperBound-designSpaceLowerBound);
labelSampleNew = design_deficit_to_label_score(designEvaluatorNew.evaluate(designSampleNew));
updatedCandidateSpace = candidateSpace.update_candidate_space(designSampleNew,labelSampleNew);
shapeDefinitionSample = updatedCandidateSpace.DesignSampleDefinition(updatedCandidateSpace.IsShapeDefinition,:);

figure;
updatedCandidateSpace.plot_candidate_space(gcf,'FaceColor','g','FaceAlpha',0.3,'EdgeColor','none');
hold all;
plot(designSample(labelSample,1),designSample(labelSample,2),'g.');
plot(designSample(~labelSample,1),designSample(~labelSample,2),'r.');
plot(designSampleNew(labelSampleNew,1),designSampleNew(labelSampleNew,2),'c.');
plot(designSampleNew(~labelSampleNew,1),designSampleNew(~labelSampleNew,2),'m.');
plot(shapeDefinitionSample(:,1),shapeDefinitionSample(:,2),'bo');
grid minor;


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

