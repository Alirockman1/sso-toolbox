function varargout = logger_objective_function(x)
    persistent baseFunction designSample objectiveValue;

    if(isempty(baseFunction))
        designSample = [];
        objectiveValue = [];
    end

    if(nargin == 0)
        varargout{1} = designSample;
        varargout{2} = objectiveValue;
        return;
    end

    if(isa(x,'function_handle'))
        baseFunction = x;
        designSample = [];
        objectiveValue = [];
        return;
    end

    if(isempty(x))
        varargout{1} = [];
        return;
    end

    if(isempty(baseFunction))
        error('No base function set');
    end
    
    [currentObjectiveValue] = baseFunction(x);
    
    designSample = [designSample; x];
    objectiveValue = [objectiveValue; currentObjectiveValue];

    varargout{1} = currentObjectiveValue;
end