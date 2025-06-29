function performanceMeasure = coil_spring_design(designSample,systemParameter)
	%% unwrap design variables
	wireDiameter = designSample(:,1); % d
	coilDiameter = designSample(:,2); % D
	nWindings = designSample(:,3); % n
	spacingWiresRelaxed = designSample(:,4); % s
	density = designSample(:,5); % rho
	shearModulus = designSample(:,6); % G
	yieldStress = designSample(:,7); % sigmaY

	%% bottom-up mappings
	% k + G*d^4/(8*D^3*n)
	stiffness = (shearModulus.*wireDiameter.^4)./(8.*coilDiameter.^3.*nWindings);
	% m = pi^2/4*rho*d^2*D*n
	mass = (pi^2/4).*density.*wireDiameter.^2*coilDiameter.*nWindings;
	% Wo = D + d
	outerDiameterOperationSpace = coilDiameter + wireDiameter;
	% Wi = D - d
	innerDiameterOperationSpace = coilDiameter - wireDiameter;
	% uc = (s-d)*n
	deformationCompaction = (spacingWiresRelaxed - wireDiameter).*nWindings;
	% uY = sigmaY*(2*pi*d^3)/(K*D) (tau_max = d/2 * MT/Jp = MT/(pi*d^3) = K*uY*D/(2*pi*d^3) =~ sigma_eq^max)
	deformationPlasticYield = yieldStress.*(2*pi.*wireDiameter.^3)./(stiffness.*coilDiameter);
	% uL = 0.1*D*n (change of angle: delta alpha =~ u/(Dn) << 0.1)
	deformationEndLinearity = 0.1.*coilDiameter.*nWindings;
	% H + L0 + uop (uop + min(uY,uL))
	lengthRelaxedState = spacingWiresRelaxed.*nWindings;
	heightOperatingSpace = lengthRelaxedState + min(deformationPlasticYield,deformationEndLinearity);


	%% wrap quantities of interest
	performanceMeasure(:,1) = stiffness; % K
	performanceMeasure(:,2) = mass; % m
	performanceMeasure(:,3) = heightOperatingSpace; % H
	performanceMeasure(:,4) = outerDiameterOperationSpace; % Wo
	performanceMeasure(:,5) = innerDiameterOperationSpace; % Wi
	performanceMeasure(:,6) = deformationCompaction; % uc
	performanceMeasure(:,7) = deformationPlasticYield; % uY
	performanceMeasure(:,8) = deformationEndLinearity; % uL
end