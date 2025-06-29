function performanceMeasure = coil_spring_design_given_material(designSample,systemParameter)	
	nSample = size(designSample,1);
	performanceMeasure = coil_spring_design([designSample,repmat(systemParameter,nSample,1)]);
end