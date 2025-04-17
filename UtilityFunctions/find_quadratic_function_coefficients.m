function [quadraticCoefficient, linearCoefficient, constantTerm] = find_quadratic_function_coefficients(functionHandle,nDimension)
%
%   Finds a,b,c assuming function can be written as f(x) = 1/2*x*a*x' + b*x' + c (where x is
%   a row vector of length nDimension)
%
%   Inputs:
%       functionHandle - function handle
%       nDimension - number of dimensions
%
%   Outputs:
%       quadraticCoefficient - quadratic coefficient
%       linearCoefficient - linear coefficient
%       constantTerm - constant term
    
    % Get constant term (function evaluated at zero)
    constantTerm = functionHandle(zeros(1,nDimension))';
    nOutput = size(constantTerm,1);

    % get diagonal and linear coefficients exploiting unit vector
    quadraticCoefficient = nan(nDimension,nDimension,nOutput);
    linearCoefficient = nan(nOutput,nDimension);
    unitVector = zeros(1,nDimension);
    for iDimension = 1:nDimension
        unitVector(iDimension) = 1;

        % 1/2*a(i,i) + b(i) = f(unitVector) - c = sum
        % 1/2*a(i,i) - b(i) = f(-unitVector) - c = diff
        % -> a(i,i) = sum + diff
        % -> b(i) = (sum - diff)/2
        halfDiagonalPlusLinearCoefficient = functionHandle(unitVector)' - constantTerm;
        halfDiagonalMinusLinearCoefficient = functionHandle(-unitVector)' - constantTerm;

        quadraticCoefficient(iDimension,iDimension,:) = halfDiagonalPlusLinearCoefficient + halfDiagonalMinusLinearCoefficient;
        linearCoefficient(:,iDimension) = (halfDiagonalPlusLinearCoefficient - halfDiagonalMinusLinearCoefficient)./2;
        
        unitVector(iDimension) = 0;
    end

    % get off-diagonal coefficients exploiting pair vector
    pairVector = zeros(1,nDimension);
    for iDimension = 1:nDimension
        pairVector(iDimension) = 1;

        for jDimension = (iDimension+1):nDimension
            pairVector(jDimension) = 1;

            % 1/2*(a(i,i)+a(i,j)+a(j,i)+a(j,j)) + b(i) + b(j) + c = f(pairVector)
            % enforcing symmetry a(i,j) = a(j,i) = f(pairVector) - b(i) - b(j) - c - 1/2*(a(i,i)+a(j,j))
            easyColumnTerm = functionHandle(pairVector)' - linearCoefficient*pairVector' - constantTerm;
            quadraticTermAsColumn = reshape(quadraticCoefficient(iDimension,iDimension,:) + quadraticCoefficient(jDimension,jDimension,:),nOutput,1);
            quadraticCoefficient(iDimension,jDimension,:) = (easyColumnTerm - 1/2*quadraticTermAsColumn);
            quadraticCoefficient(jDimension,iDimension,:) = quadraticCoefficient(iDimension,jDimension,:);
            
            pairVector(jDimension) = 0; 
        end

        pairVector(iDimension) = 0;
    end
end
