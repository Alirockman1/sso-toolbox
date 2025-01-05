import numpy

def car_crash_2d_python(designSample,systemParameter):
    """
    Evaluate a simplified 2D car crash model and return key performance 
    measures.

    This function calculates three quantities for each design sample:
    
    1. Energy Remaining: The energy left in the vehicle after deformation
    2. Maximum Acceleration: The maximum acceleration the vehicle experiences.
    3. Deformation Order: A measure of which component deforms first
       
    The function returns these three performance measures as a NumPy array of 
    shape (n_samples, 3), where each row corresponds to one design sample.

    Parameters
    ----------
    designSample : array_like
        A 2D array or list of shape (nSample, 2) specifying design samples. 
        Each row is one design sample consisting of:
        - force1 (F1): The collapse load for component 1 (column 0).
        - force2 (F2): The collapse load for component 2 (column 1).
    
    systemParameter : array_like
        A 1D array or list of shape (4,) specifying parameters of the system:
        - vehicleMass (m) :float
            Vehicle mass.
        - criticalDisplacement1 (d1c) :float
            Critical displacement for component 1.
        - criticalDisplacement2 (d2c) :float
            Critical displacement for component 2.
        - initialVehicleSpeed (v0) :float
            Initial vehicle speed.

    Returns
    -------
    performanceMeasure : numpy.ndarray
        A 2D NumPy array of shape (n_samples, 3) where each row contains:
        - energyRemaining :float
            Remaining energy in the vehicle after deformation.
        - maximumAcceleration :float
            The maximum acceleration experienced by the vehicle.
        - deformationOrder :float
            A measure indicating which component deforms first 
            (positive if F1 > F2, negative if F2 > F1).

    Notes
    -----
    - The function internally converts inputs to NumPy arrays.
    - All forces, mass, speeds, and displacements are assumed to be given in 
      compatible units such that the energy computations are valid (e.g., 
      Newtons for force, meters for displacement, kg for mass, m/s for speed).
    """

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