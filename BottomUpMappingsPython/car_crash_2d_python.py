import numpy

def car_crash_2d_python(designSample,systemParameter):
    designSample = numpy.array(designSample,ndmin=2)
    systemParameter = numpy.array(systemParameter)

    # Unwrap inputs
    force1 = designSample[:, 0] # F1 
    force2 = designSample[:, 1] # F2 
    vehicleMass = systemParameter[0] # m 
    criticalDisplacement1 = systemParameter[1] # d1c 
    criticalDisplacement2 = systemParameter[2] # d2c 
    initialVehicleSpeed = systemParameter[3] # v0 
    
    # Energy remaining = energy vehcile has - deformation energy
    # Erem = 1/2*m*v0^2 - (F1*d1c + F2*d2c)
    energyRemaining = (1/2)*vehicleMass*(initialVehicleSpeed**2) - (force1*criticalDisplacement1 + force2*criticalDisplacement2)
    
    # Maximum Acceleration = force of strongest collapse load / total mass
    # amax = F2/m
    maximumAcceleration = force2/vehicleMass
    
    # Deformation order = force of component 1 (lesser) - force of component 2 (greater)
    # order = F1 - F2
    deformationOrder = force1 - force2
    
    # Wrap outputs
    performanceMeasure = numpy.transpose(numpy.vstack((energyRemaining,maximumAcceleration,deformationOrder)))
    physicalFeasibilityMeasure = numpy.zeros(performanceMeasure.shape[0])
    
    return performanceMeasure
    # return (performanceMeasure,physicalFeasibilityMeasure)