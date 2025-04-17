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

    currentEqualityConstraintValue = [];
    if(nargout(baseFunction) == 1)
        currentInequalityConstraintValue = baseFunction(x);
    elseif(nargout(baseFunction) == 2)
        [currentInequalityConstraintValue,currentEqualityConstraintValue] = baseFunction(x);
    else
        try
            [currentInequalityConstraintValue,currentEqualityConstraintValue] = baseFunction(x);
        catch
            try
                currentInequalityConstraintValue = baseFunction(x);
            catch
                error('logger_constraint_function:InvalidNumberOfOutputs', 'The function must return either one or two outputs.');
            end
        end
    end

    designSample = [designSample; x];
    inequalityConstraintValue = [inequalityConstraintValue; currentInequalityConstraintValue];
    equalityConstraintValue = [equalityConstraintValue; currentEqualityConstraintValue];

    varargout{1} = currentInequalityConstraintValue;
    varargout{2} = currentEqualityConstraintValue;
end