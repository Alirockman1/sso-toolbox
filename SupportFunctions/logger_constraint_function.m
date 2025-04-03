function varargout = logger_constraint_function(x)
    persistent baseFunction designSample inequalityConstraintValue equalityConstraintValue;

    if(isempty(baseFunction))
        designSample = [];
        inequalityConstraintValue = [];
        equalityConstraintValue = [];
    end

    if(nargin == 0)
        varargout{1} = designSample;
        varargout{2} = inequalityConstraintValue;
        varargout{3} = equalityConstraintValue;
        return;
    end

    if(isa(x,'function_handle'))
        baseFunction = x;
        designSample = [];
        inequalityConstraintValue = [];
        equalityConstraintValue = [];
        return;
    end

    if(isempty(x))
        varargout{1} = [];
        return;
    end

    if(isempty(baseFunction))
        error('No base function set');
    end

    if(nargout(baseFunction) == 1)
        currentInequalityConstraintValue = baseFunction(x);
        currentEqualityConstraintValue = [];
    else %if(nargout(baseFunction) == 2)
        [currentInequalityConstraintValue, currentEqualityConstraintValue] = baseFunction(x);
    end
    designSample = [designSample; x];
    inequalityConstraintValue = [inequalityConstraintValue; currentInequalityConstraintValue];
    equalityConstraintValue = [equalityConstraintValue; currentEqualityConstraintValue];

    varargout{1} = currentInequalityConstraintValue;
    varargout{2} = currentEqualityConstraintValue;
end