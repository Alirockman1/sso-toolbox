function [plotData,evaluationData] = evaluate_selective_design_space_projection(designEvaluator,designBox,designSpaceLowerBound,designSpaceUpperBound,desiredPairs,varargin)

    %% Parse Inputs
    p = inputParser;
    addOptional(p,'SamplingMethod',@sampling_latin_hypercube);
    addOptional(p,'SamplingOptions',{},@iscell);
    addOptional(p,'NumberSamplesPerPlot',3000,@isnumeric);
    addOptional(p,'MarkerColorsCriterion','random',@(x) ismember(x,{'random','worst'}));
    parse(p,varargin{:});
    options = p.Results;
    
    % Gather all evaluations necessary for all the plots
    nPlot = size(desiredPairs,1);
    nDesignVariable = size(designSpaceLowerBound,2);
    nPointTotal = nPlot*options.NumberSamplesPerPlot;
    designSample = nan(nPointTotal,nDesignVariable);

    %% Random Sampling for Relevant Scatter Plot
    for j=1:nPlot
        % Calculate row indices for this plot's samples
        startRow = (j-1)*options.NumberSamplesPerPlot + 1;
        endRow = j*options.NumberSamplesPerPlot;
        currentRow = startRow:endRow;
        
        currentAxes = [desiredPairs(j,1),desiredPairs(j,2)];
        otherAxes = ~ismember(1:nDesignVariable,currentAxes);
        
        scatterBound = nan(2, nDesignVariable);
        scatterBound(:,currentAxes) = [designSpaceLowerBound(currentAxes);designSpaceUpperBound(currentAxes)];
        scatterBound(:,otherAxes) = designBox(:,otherAxes);

        designSample(currentRow,:) = options.SamplingMethod(scatterBound, options.NumberSamplesPerPlot, options.SamplingOptions{:});
    end


    %% Evaluate Samples
    [performanceDeficit,physicalFeasibilityDeficit,evaluatorOutput] = designEvaluator.evaluate(designSample);

    if(isempty(performanceDeficit))
        isGoodPerformance = true(size(designSample,1),1);
    else
        isGoodPerformance = design_deficit_to_label_score(performanceDeficit);
    end

    if(isempty(physicalFeasibilityDeficit))
        isPhysicallyFeasible = true(size(designSample,1),1);
    else
        isPhysicallyFeasible = design_deficit_to_label_score(physicalFeasibilityDeficit);
    end

    %% Tag Bad Designs w.r.t. which requirement it violates
    iRequirementViolated = zeros(nPointTotal,1);
    if(strcmpi(options.MarkerColorsCriterion,'random'))
        for i=1:nPointTotal
            if(~isGoodPerformance(i))
                currentRequirementViolated = find(performanceDeficit(i,:)>0);
                randomizedViolatedRequirement = currentRequirementViolated(randperm(length(currentRequirementViolated)));
                iRequirementViolated(i) = randomizedViolatedRequirement(1);
            end
        end
    else % 'worst'
        [~,iWorstDeficit] = max(performanceDeficit(~isGoodPerformance,:),[],2);
        iRequirementViolated(~isGoodPerformance) = iWorstDeficit;
    end

    plotData = struct(...
        'DesignSample',[],...
        'RequirementViolated',[],...
        'IsPhysicallyFeasible',[],...
        'DesignSpaceLowerBound',[],...
        'DesignSpaceUpperBound',[],...
        'DesignBox',[]);
    
    for j=1:nPlot
        % Calculate row indices for this plot's samples
        startRow = (j-1)*options.NumberSamplesPerPlot + 1;
        endRow = j*options.NumberSamplesPerPlot;
        currentRow = startRow:endRow;

        currentPair = [desiredPairs(j,1),desiredPairs(j,2)];

        % Store data for this plot
        plotData(j) = struct(...
            'DesignSample',designSample(currentRow,currentPair),...
            'RequirementViolated',iRequirementViolated(currentRow),...
            'IsPhysicallyFeasible',isPhysicallyFeasible(currentRow),...
            'DesignSpaceLowerBound',designSpaceLowerBound(currentPair),...
            'DesignSpaceUpperBound',designSpaceUpperBound(currentPair),...
            'DesignBox',designBox(:,currentPair));
    end

    if(nargout>1)
        evaluationData.DesignSample = designSample;
        evaluationData.PerformanceDeficit = performanceDeficit;
        evaluationData.PhysicalFeasibilityDeficit = physicalFeasibilityDeficit;
        evaluationData.EvaluatorOutput = evaluatorOutput;
    end
end