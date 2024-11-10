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


%% grow candidate space
grownCandidateSpace = candidateSpace.grow_candidate_space(0.1);

figure;
plot(designSample(labelSample,1),designSample(labelSample,2),'g.');
hold all;
plot(designSample(~labelSample,1),designSample(~labelSample,2),'r.');
grownCandidateSpace.plot_candidate_space(gcf,'FaceColor','g','FaceAlpha',0.5,'EdgeColor','none');
grid minor;


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

