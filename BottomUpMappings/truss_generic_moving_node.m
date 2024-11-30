function performanceMeasure = truss_generic_moving_node(designSample,systemParameter)


	% check which entries are design variables based on nan entry
	isDesignVariable = isnan(systemParameter.BaseNodePosition);
	isTrussTip = (systemParameter.NodeForce~=0);
	nDimension = size(systemParameter.BaseNodePosition,2);
	nElement = size(systemParameter.NodeElement,1);

	nodePosition = systemParameter.BaseNodePosition;
	nSample = size(designSample,1);
	performanceMeasure = nan(nSample,2+nElement);
	for i=1:nSample
		nodePosition(isDesignVariable) = column_vector_to_row_major_matrix(designSample(i,:)',nDimension);

		[nodeDisplacement,~,elementAxialForce] = ...
			truss_analysis(...
				nodePosition,...
				systemParameter.FixedDegreesOfFreedom,...
				systemParameter.NodeForce,...
				systemParameter.NodeElement,...
				systemParameter.ElementCrossSectionArea,...
				systemParameter.ElementYoungsModulus);
		elementStress = truss_deformed_stress(elementAxialForce,systemParameter.ElementCrossSectionArea);

		nodeDistance = nodePosition(systemParameter.NodeElement(:,2),:)-nodePosition(systemParameter.NodeElement(:,1),:);
	    elementLength = vecnorm(nodeDistance,2,2);
		totalMass = sum(systemParameter.ElementDensity.*elementLength.*systemParameter.ElementCrossSectionArea);

        % deal with numerical errors
        displacement = -nodeDisplacement(isTrussTip);
        displacement(isnan(displacement)) = +inf;
        elementStress(isnan(elementStress)) = +inf;

		performanceMeasure(i,:) = [displacement,totalMass,abs(elementStress')];
	end
end

