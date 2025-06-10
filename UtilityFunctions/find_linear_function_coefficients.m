function [linearCoefficient, constantTerm] = find_linear_function_coefficients(functionHandle,nDimension)
%
%   Finds a,b assuming function can be written as f(x) = a*x' + b (where x is
%   a row vector of length nDimension)
%
%   Inputs:
%       functionHandle - function handle
%       nDimension - number of dimensions
%
%   Outputs:
%       linearCoefficient - linear coefficient
%       constantTerm - constant term

    % Get constant term (function evaluated at zero)
    constantTerm = functionHandle(zeros(1,nDimension))';
    nOutput = size(constantTerm,1);

    % Get linear coefficients (function evaluated at unit vector for each dimension minus constant term)
    linearCoefficient = nan(nOutput,nDimension);
    unitVector = zeros(1,nDimension);
    for iDimension = 1:nDimension
        unitVector(iDimension) = 1;

        linearCoefficient(:,iDimension) = functionHandle(unitVector)' - constantTerm;
        
        unitVector(iDimension) = 0;
    end
end