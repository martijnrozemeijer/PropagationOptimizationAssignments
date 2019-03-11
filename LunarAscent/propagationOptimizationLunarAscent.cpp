/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h>
#include <Tudat/Mathematics/GeometricShapes/capsule.h>
#include <Tudat/Mathematics/Statistics/randomVariableGenerator.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "../applicationOutput.h"

using namespace tudat;
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::aerodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::mathematical_constants;
using namespace tudat::reference_frames;
using namespace tudat::estimatable_parameters;

/*!
 *  Class to compute the thrust direction and magnitude for the lunar ascent vehicle. The current inputs set a
 *  constant thrust magnitude, and a thrust direction in y-(-z) plane of the vertical frame define linearly in time using
 *  equispaced nodes. These settings are to be modified for the assignment.
 */
class LunarAscentThrustGuidance
{
public:

    //! Contructor
    /*!
     * Contructor
     * \param vehicleBody Body object for the ascent vehicle
     * \param initialTime Start time of the propagation
     * \param parameterVector Vector of independent variables to be used for thrust parameterization:
     *   - Entry 0: Constant thrust magnitude
     *   - Entry 1: Constant spacing in time between nodes
     *   - Entry 2-6: Thrust angle theta, at five nodes
     */
    LunarAscentThrustGuidance(
            const std::shared_ptr< Body > vehicleBody,
            const double initialTime,
            const std::vector< double > parameterVector ):
        vehicleBody_( vehicleBody ),
        parameterVector_( parameterVector )
    {
        // Retrieve parameters of thrust profile
        thrustMagnitude_ = parameterVector_.at( 0 );
        timeInterval_ = parameterVector_.at( 1 );

        // Create interpolator for thrust angle
        double currentTime = initialTime;
        for( unsigned int i = 0; i < parameterVector_.size( ) - 2; i++ )
        {
            thrustAngleMap_[ currentTime ] = parameterVector_.at( i + 2 );
            currentTime += timeInterval_;
        }
        thrustAngleInterpolator_ = createOneDimensionalInterpolator(
                    thrustAngleMap_, std::make_shared< InterpolatorSettings >(
                        linear_interpolator, huntingAlgorithm, false, use_boundary_value ) );
    }

    //! Function that computes the inertial thrust direction for each state derivative function evaluation
    Eigen::Vector3d getCurrentThrustDirection( const double currentTime )
    {
        // Retrieve thrust angle
        double currentThrustAngle = thrustAngleInterpolator_->interpolate( currentTime );

        // Set thrust in V-frame
        Eigen::Vector3d thrustDirectionInVerticalFrame =
                ( Eigen::Vector3d( ) << 0.0, std::sin( currentThrustAngle ), -std::cos( currentThrustAngle ) ).finished( );

        // Retrieve rotation from V-frame to inertial frame
        Eigen::Quaterniond verticalToInertialFrame =
                vehicleBody_->getFlightConditions( )->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                    vertical_frame, inertial_frame );

        // Return thrust direction
        return verticalToInertialFrame * thrustDirectionInVerticalFrame;

    }

    //! Function that computes the thrust magnitude for each state derivative function evaluation
    double getCurrentThrustMagnitude( const double currentTime )
    {
        return thrustMagnitude_;
    }

private:

    //! Object containing properties of the vehicle
    std::shared_ptr< Body > vehicleBody_;

    //! Parameter containing the solution parameter vector
    std::vector< double > parameterVector_;

    //! Map containing the thrust (value) as a function of time (key)
    std::map< double, double > thrustAngleMap_;

    //! Object that interpolates the thrust as a function of time
    std::shared_ptr< OneDimensionalInterpolator< double, double > > thrustAngleInterpolator_;

    //! Constant time between thrust angle nodes
    double timeInterval_;

    //! Constant magnitude of the thrust
    double thrustMagnitude_;

};

/*!
 *   This function computes the dynamics of a lunar ascent vehicle, starting at zero velocity on the Moon's surface. The only
 *   accelerations acting on the spacecraft are the Moon's point-mass gravity, and the thrust of the vehicle. Both the
 *   translational dynamics and mass of the vehicle are propagated, using a fixed specific impulse.
 *
 *   The thrust is computed based on a fixed thrust magnitude, and a variable thrust direction. The trust direction is defined
 *   on a set of 5 nodes, spread eavenly in time. At each node, a thrust angle theta is defined, which gives the angle between the
 *   -z and y angles in the ascent vehicle's vertical frame (see Mooij, 1994). Between the nodes, the thrust is linearly
 *   interpolated. If the propagation goes beyond the bounds of the nodes, the boundary value is used.
 *
 *   The propagation is terminated as soon as one of the following conditions is met:
 *
 *   - Altitude > 100 km
 *   - Altitude < 0 km
 *   - Propagation time > 3600 s
 *   - Vehicle mass < 1000 kg
 *
 *   Key outputs:
 *
 *   propagatedStateHistory Numerically propagated Cartesian state
 *   dependentVariableHistory Dependent variables saved during the state propagation of the ascent *
 *
 *   Input parameters:
 *
 *   thrustParameters: Vector contains the following:
 *
 *   - Entry 0: Constant thrust magnitude
 *   - Entry 1: Constant spacing in time between nodes
 *   - Entry 2-6: Thrust angle theta, at five nodes
 */
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

//    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_pa_de421_1900-2050.bpc" );
//    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_080317.tf" );
//    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_assoc_pa.tf" );


    std::string outputPath = tudat_applications::getOutputPath( "LunarAscent" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define vehicle settings
    double vehicleMass = 4.7E3;
    double vehicleDryMass = 2.25E3;
    double constantSpecificImpulse = 311.0;

    // Define simulation settings
    double initialTime = 0.0;
    double maximumDuration = 86400.0;
    double terminationAltitude = 100.0E3;

    // Define initial spherical elements for vehicle.
    Eigen::Vector6d ascentVehicleSphericalEntryState;
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Moon" ) + 100.0;
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 0.6875 );
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 23.4333 );
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 0.0;
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            unit_conversions::convertDegreesToRadians( 90.0 );
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) =
            unit_conversions::convertDegreesToRadians( 90.0 );

    // Convert vehicle state from spherical elements to body-fixed Cartesian elements.
    Eigen::Vector6d bodyFixedSystemInitialState = convertSphericalOrbitalToCartesianState(
                ascentVehicleSphericalEntryState );

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > thrustParameters =
    {15727.8404745739,74.82920403592288,-0.09481475441716612,0.2081478221807629,0.06612515454180534,0.7414646069519222,-0.4368525084573776};

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT                   //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create solar system bodies
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Earth" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
                    getDefaultBodySettings( bodiesToCreate, -600.0, 950.0 );

//    std::string moonFrame;
//    for (unsigned int rotationModelCase = 0; rotationModelCase <2; rotationModelCase++ )
//    if( rotationModelCase == 0 )
//    {
//        moonFrame = "IAU_Moon_Simplified";
//    }
//    else if( rotationModelCase == 1 )
//    {
//        moonFrame = "IAU_Moon";
//    }
//    if( rotationModelCase == 2 )
//    {
//        moonFrame = "MOON_ME";
//    }
//    else if( rotationModelCase == 2 )
//    {
//        moonFrame = "MOON_PA";
//    }

//    if( rotationModelCase == 0 )
//    {
//        bodySettings[ "Moon" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
//                    "ECLIPJ2000", moonFrame,
//                    spice_interface::computeRotationQuaternionBetweenFrames(
//                        "ECLIPJ2000", "IAU_Moon", initialTime ),
//                        initialTime, spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame(
//                            "ECLIPJ2000", "IAU_Moon", initialTime ).norm( ) );
//    }
//    else
//    {

//        bodySettings[ "Moon" ]->rotationModelSettings = std::make_shared< RotationModelSettings >(
//                    spice_rotation_model, "ECLIPJ2000", moonFrame );

//    }
//    bodySettings[ "Earth" ]->rotationModelSettings = std::make_shared< RotationModelSettings >(
//                            spice_rotation_model, "ECLIPJ2000", moonFrame );

//    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrientation( "ECLIPJ2000" );
//    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "ECLIPJ2000" );

//    std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
//                bodySettings[ "Moon" ]->gravityFieldSettings )->resetAssociatedReferenceFrame( moonFrame );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

//    double referenceAreaRadiation = 15.0;
//    double radiationPressureCoefficient = 0.01;
//    std::vector< std::string > occultingBodies;
//    occultingBodies.push_back( "Moon" );
//    std::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings =
//            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
//                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
//    bodyMap[ "Vehicle" ]->setRadiationPressureInterface(
//                "Sun", createRadiationPressureInterface(
//                    vehicleRadiationPressureSettings, "Vehicle", bodyMap ) );

    // Finalize body creation.


    setGlobalFrameBodyEphemerides( bodyMap, "Moon", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define thrust functions
    std::shared_ptr< LunarAscentThrustGuidance > thrustGuidance =
            std::make_shared< LunarAscentThrustGuidance >(
                bodyMap.at( "Vehicle" ), initialTime, thrustParameters );
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
            std::bind( &LunarAscentThrustGuidance::getCurrentThrustDirection, thrustGuidance, std::placeholders::_1 );
    std::function< double( const double ) > thrustMagnitudeFunction =
            std::bind( &LunarAscentThrustGuidance::getCurrentThrustMagnitude, thrustGuidance, std::placeholders::_1 );

    std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< CustomThrustDirectionSettings >( thrustDirectionFunction );
    std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, [ = ]( const double ){ return constantSpecificImpulse; } );

    // Vary acceleration settings
//    for (int i = 1; i <8; i++){
    // Define acceleration settings
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                       thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 )  );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );

//    DECOMMENT FOR ASSIGNMENT 2,1
//    //Create points mass for first propagation
//    if (i >0 && i < 4){
//        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
//                                                        basic_astrodynamics::central_gravity ) );
//    }
//    //Include Earth as point mass
//    if (i > 1 && i < 5){
//        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
//                                                    basic_astrodynamics::central_gravity ) );
//     }
//    //Include Sun as point mass
//    if ( i > 2 ){
//        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
//                                                    basic_astrodynamics::central_gravity ) );
//    }
//    //Include spherical harmonics for Moon
//    if (i == 4 ){
//        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 )  );
//     }
//    //Include spherical harmonics for Moon
//    if (i > 4 ){
//        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 )  );
//     }

//    //Include cannon radiation pressure from Sun
//    if (i > 5){
//        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
//                                                       basic_astrodynamics::cannon_ball_radiation_pressure ) );
//    }
//    if (i > 6){
//        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 )  );

//    }
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Moon" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Convert the state to the global (inertial) frame.
    std::shared_ptr< ephemerides::RotationalEphemeris > moonRotationalEphemeris =
            bodyMap.at( "Moon" )->getRotationalEphemeris( );
    Eigen::VectorXd systemInitialState = transformStateToGlobalFrame(
                bodyFixedSystemInitialState, initialTime, moonRotationalEphemeris );

    // Create termination settings.
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           initialTime + maximumDuration ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               altitude_dependent_variable, "Vehicle", "Moon" ), terminationAltitude, false ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               altitude_dependent_variable, "Vehicle", "Moon" ), 0.0, true ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               current_body_mass_dependent_variable, "Vehicle" ), vehicleDryMass, true ) );
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = std::make_shared<
            PropagationHybridTerminationSettings >( terminationSettingsList, true );

    // Define dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Vehicle", "Moon" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_speed_dependent_variable, "Vehicle", "Moon" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Vehicle", flight_path_angle ) );

    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );
    TranslationalPropagatorType propagatorType;
    propagatorType = cowell;

    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = (
                createMassRateModel( "Vehicle", std::make_shared< FromThrustMassModelSettings >( 1 ),
                                     bodyMap, accelerationModelMap ) );
    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Vehicle" }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );
    // Define translational state propagation settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                terminationSettings, propagatorType );

    // Define full propagation settings
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector =
    { translationalStatePropagatorSettings, massPropagatorSettings };
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

    // Define time to compare to constant step size 5.0 for propagation/integration schemes
//    std::cout << "Use RK4 Integrator to get 5 second step size"<< std::endl;
//    std::shared_ptr< IntegratorSettings< > > integratorSettingsComparitor;
//    double fixedStepSize = 5.0;
//    integratorSettingsComparitor =
//        std::make_shared< IntegratorSettings < > >
//            ( rungeKutta4, initialTime, fixedStepSize );
//    SingleArcDynamicsSimulator< > dynamicsSimulatorComparitor(
//                bodyMap, integratorSettingsComparitor, propagatorSettings );
//    std::map< double, Eigen::VectorXd > propagatedStateHistoryComparitor = dynamicsSimulatorComparitor.getEquationsOfMotionNumericalSolution( );

    std::shared_ptr< IntegratorSettings< > > integratorSettings;
    double relativeErrorToleranceRK87 = 1e-12;
    double absoluteErrorToleranceRK87 = 1e-12;
    double initialTimeStep = 1.0;
    double minimumStepSize = 0.000001;
    double maximumStepSize = 5.0;
    integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
              rungeKuttaVariableStepSize , initialTime,  initialTimeStep
              , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKutta87DormandPrince
              , minimumStepSize, maximumStepSize, absoluteErrorToleranceRK87, relativeErrorToleranceRK87);

//    SingleArcDynamicsSimulator< > dynamicsSimulator(
//                bodyMap, integratorSettings, propagatorSettings );
//    std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
//    std::map< double, unsigned int > integratorData = dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Vehicle", systemInitialState, "Moon" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter  ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 2, "Moon", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 2, "Moon", spherical_harmonics_sine_coefficient_block ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap, accelerationModelMap);

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
                std::shared_ptr< IntegratorSettings< double > >( ) , false, true );

    std::map< double, Eigen::VectorXd > nominalIntegrationResult =
            variationalEquationsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::MatrixXd > stateTransitionResult =
            variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 0 );
    std::map< double, Eigen::MatrixXd > sensitivityResult =
            variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 1 );

//    // Define interpolator settings
//    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
//            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 );
//    // Create interpolation of integration/propagator
//    std::map< double, Eigen::VectorXd > InterpolatedResult;
//    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > comparitorInterpolator =
//            interpolators::createOneDimensionalInterpolator(
//                propagatedStateHistory, interpolatorSettings );
//    for( auto stateIterator : propagatedStateHistoryComparitor )
//    {
//        double currentTime = stateIterator.first;
//        InterpolatedResult[ currentTime ] =
//                comparitorInterpolator->interpolate( currentTime );
//    }
    input_output::writeDataMapToTextFile( nominalIntegrationResult, "integrationresult.dat", outputPath );
    input_output::writeDataMapToTextFile( stateTransitionResult, "stateTransitionResult.dat", outputPath );
    input_output::writeDataMapToTextFile( sensitivityResult, "sensitivityResult.dat", outputPath );

//    input_output::writeDataMapToTextFile( dependentVariableHistory, "dependentVariables.dat", outputPath );
//    input_output::writeDataMapToTextFile( integratorData, "integratorData.dat", outputPath );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             ANALYSE UNCERTANTIES IN MODELS/STATE              /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//        // Integrator settings
//        std::shared_ptr< IntegratorSettings< > > integratorSettings;
//        double initialTimeStep = 1.0;
//        double minimumStepSize = 0.000001;
//        double maximumStepSize = 5;
//        std::string behindText;
//        for (int j = 1; j<3; j++){
//            if (j == 1){
//                double relativeErrorToleranceRK45 = 1e-12;
//                double absoluteErrorToleranceRK45 = 1e-12;
//                behindText = "-12";
//                integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                          rungeKuttaVariableStepSize , initialTime,  initialTimeStep
//                          , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg45
//                          , minimumStepSize, maximumStepSize, relativeErrorToleranceRK45, absoluteErrorToleranceRK45);

//            }
//            else if (j == 2){
//                double relativeErrorToleranceRK87 = 1e-12;
//                double absoluteErrorToleranceRK87 = 1e-12;
//                behindText = "-12";
//                integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                          rungeKuttaVariableStepSize , initialTime,  initialTimeStep
//                          , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKutta87DormandPrince
//                          , minimumStepSize, maximumStepSize, absoluteErrorToleranceRK87, relativeErrorToleranceRK87);

//            }
//            // Compute dynamics of integrator/propagator combinations
//            SingleArcDynamicsSimulator< > dynamicsSimulator(
//                        bodyMap, integratorSettings, propagatorSettings );
//            std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//            std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
//            std::map< double, unsigned int > integratorData = dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations();


//            // Define interpolator settings
//            std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
//                    std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 );
//            // Create interpolation of integration/propagator
//            std::map< double, Eigen::VectorXd > InterpolatedResult;
//            std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > comparitorInterpolator =
//                    interpolators::createOneDimensionalInterpolator(
//                        propagatedStateHistory, interpolatorSettings );
//            for( auto stateIterator : propagatedStateHistoryComparitor )
//            {
//                double currentTime = stateIterator.first;
//                InterpolatedResult[ currentTime ] =
//                        comparitorInterpolator->interpolate( currentTime );
//            }
//        input_output::writeDataMapToTextFile( InterpolatedResult, "interpolatedStateHistory" + std::to_string(i) + std::to_string(j) + "_" + behindText+ ".dat", outputPath );
//        input_output::writeDataMapToTextFile( integratorData, "integratorData" + std::to_string(i) + std::to_string(j) + "_" + behindText+ ".dat", outputPath );

//        }
//    }


//        else if (i == 3){
//        std::cout<< "gauss propagation" << std::endl;
//        propagatorType = gauss_keplerian;
//        }
//        else if (i == 4){
//        std::cout<< "quaternions propagation" << std::endl;
//        propagatorType = unified_state_model_quaternions;
//        }
//        else if (i == 5){
//        std::cout<< "rodrigues parameters propagation" << std::endl;
//        propagatorType = unified_state_model_modified_rodrigues_parameters;
//        }
//        else if (i == 6){
//        std::cout<< "exponential_map propagation" << std::endl;
//        propagatorType = unified_state_model_exponential_map;
//        }
    // KEEP FOR ASSIGNMENT 1



////     STOP KEEPING FOR ASSIGNMENT 2
//    std::cout << "RK4 integration method"<< std::endl;
//    std::shared_ptr< IntegratorSettings< > > integratorSettings;
//    double fixedStepSize = 1.0;
//    integratorSettings =
//        std::make_shared< IntegratorSettings< > >
//            ( rungeKutta4, initialTime, fixedStepSize );


//        input_output::writeDataMapToTextFile( interpolatedBenchmarkResult, "interpolatedStateHistory" + std::to_string(i) + "_" + behindText+ ".dat", outputPath );
//        input_output::writeDataMapToTextFile( integratorData, "integratorData" + std::to_string(i) + "_" + behindText + ".dat", outputPath );



        // HERE STARTS ASSIGNMENT 1 UNCOMMENT FOR USE --> FOR LOOP IS IMPLEMENTED
//        // Define interpolator settings
//        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
//                std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 );

//    // Create benchmark results with RK78 10^-14 - 1 second initial timestep
//    std::cout << "Creating benchmark data RK78"<< std::endl;
//    std::shared_ptr< IntegratorSettings< > > integratorSettingsBenchMarkRK78;
//    integratorSettingsBenchMarkRK78 = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//              rungeKuttaVariableStepSize , initialTime,  1
//              , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78
//              , 0.000001, 5, 1e-12, 1e-12);
//    SingleArcDynamicsSimulator< > dynamicsSimulatorRK78(
//                bodyMap, integratorSettingsBenchMarkRK78, propagatorSettings );
//    std::map< double, Eigen::VectorXd > propagatedHistoryBenchMark78 = dynamicsSimulatorRK78.getEquationsOfMotionNumericalSolution( );
//    input_output::writeDataMapToTextFile( propagatedHistoryBenchMark78, "stateHistoryBenchMark78.dat", outputPath );

    // Benchmark creation, selected through comparison
//    std::cout << "Creating benchmark data RK87"<< std::endl;
//    std::shared_ptr< IntegratorSettings< > > integratorSettingsBenchMarkRK87;
//    integratorSettingsBenchMarkRK87 = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//              rungeKuttaVariableStepSize , initialTime,  1
//              , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKutta87DormandPrince
//              , 0.000001, 5, 1e-12, 1e-12);
//    SingleArcDynamicsSimulator< > dynamicsSimulatorRK87(
//                bodyMap, integratorSettingsBenchMarkRK87, propagatorSettings );
//    std::map< double, Eigen::VectorXd > propagatedHistoryBenchMark = dynamicsSimulatorRK87.getEquationsOfMotionNumericalSolution( );
//    input_output::writeDataMapToTextFile( propagatedHistoryBenchMark, "stateHistoryBenchMark87.dat", outputPath );

//        // Define interpolator settings
//        std::map< double, Eigen::VectorXd > interpolatedBenchmarkResult;
//        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > benchmarkInterpolator =
//                interpolators::createOneDimensionalInterpolator(
//                    propagatedHistoryBenchMark78, interpolatorSettings );

//        for( auto stateIterator : propagatedHistoryBenchMark87 )
//        {
//            double currentTime = stateIterator.first;
//            interpolatedBenchmarkResult[ currentTime ] =
//                    benchmarkInterpolator->interpolate( currentTime );
//        }

//        input_output::writeDataMapToTextFile( interpolatedBenchmarkResult, "interpolatedPropagatedHistoryBenchMark78.dat", outputPath );

//    // Define integration settings in for loop (8 integration methods)
//    // For all
//    const double fixedStepSize = 10.0;
//    double initialTimeStep = 1.0;
//    double minimumStepSize = 0.000001;
//    double maximumStepSize = 5;
//    double relativeErrorTolerance = 1e-12;
//    double absoluteErrorTolerance = 1e-12;
//    std::string behindText;
//    // For ADM
//    double minimumOrder = 6;
//    double maximumOrder = 11;
//    std::shared_ptr< IntegratorSettings< > > integratorSettings;
//    for (int i = 5; i < 7; i++) {
//        behindText = "10";
//        if (i == 1){
//            std::cout << "Euler integration method"<< std::endl;
//            integratorSettings =
//                std::make_shared< IntegratorSettings< > >
//                    ( euler, initialTime, fixedStepSize );
//            behindText = std::to_string(fixedStepSize);
//       }
//        else if (i == 2){
//            std::cout << "RK4 integration method"<< std::endl;
//            integratorSettings =
//                std::make_shared< IntegratorSettings< > >
//                    ( rungeKutta4, initialTime, fixedStepSize );
//            behindText = std::to_string(fixedStepSize);
//        }
//        else if (i > 2 && i < 7){
//            if (i == 3){
//                std::cout << "RK45 integration method"<< std::endl;
//                double relativeErrorToleranceRK45 = 1e-12;
//                double absoluteErrorToleranceRK45 = 1e-12;
//                behindText = "12";

//                integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                          rungeKuttaVariableStepSize , initialTime,  initialTimeStep
//                          , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg45
//                          , minimumStepSize, maximumStepSize, relativeErrorToleranceRK45, absoluteErrorToleranceRK45);
//        }
//            if (i == 4){
//                std::cout << "RK56 integration method"<< std::endl;
//                integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                          rungeKuttaVariableStepSize , initialTime,  initialTimeStep
//                          , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg56
//                          , minimumStepSize, maximumStepSize, relativeErrorTolerance, absoluteErrorTolerance);
//        }
//            if (i == 5){
//                std::cout << "RK78 integration method"<< std::endl;
//                integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                          rungeKuttaVariableStepSize , initialTime,  initialTimeStep
//                          , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78
//                          , minimumStepSize, maximumStepSize, relativeErrorTolerance, absoluteErrorTolerance);
//        }
//            if (i == 6) {
//                std::cout << "RK8(7) integration method"<< std::endl;
//                integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                          rungeKuttaVariableStepSize , initialTime,  initialTimeStep
//                          , numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKutta87DormandPrince
//                          , minimumStepSize, maximumStepSize, relativeErrorTolerance, absoluteErrorTolerance);
//        }
//        }

//        else if (i == 7){
//            std::cout << "BulirschStoer Extrapolation method"<< std::endl;
//            integratorSettings =
//                std::make_shared< BulirschStoerIntegratorSettings< > >
//                    (initialTime,
//                     initialTimeStep,
//                     bulirsch_stoer_sequence,
//                     6,
//                     minimumStepSize,
//                     maximumStepSize,
//                     relativeErrorTolerance,
//                     absoluteErrorTolerance );
//        }
//        else if (i == 8){
//         std::cout << "ADM integration method"<< std::endl;
//         double relativeErrorToleranceADM = 1e-8;
//         double absoluteErrorToleranceADM = 1e-8;
//         behindText = "8";
//         minimumStepSize = 0.000001;
//         integratorSettings =
//             std::make_shared< AdamsBashforthMoultonSettings< > >
//                     ( initialTime,
//                       initialTimeStep,
//                       minimumStepSize,
//                       maximumStepSize,
//                       relativeErrorToleranceADM,
//                       absoluteErrorToleranceADM,
//                       minimumOrder,
//                       maximumOrder );

//        }
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//        // Create simulation object and propagate dynamics.
//        SingleArcDynamicsSimulator< > dynamicsSimulator(
//                    bodyMap, integratorSettings, propagatorSettings );

//        std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//        std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
//        std::map< double, unsigned int > integratorData = dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations();

//        // Define interpolator settings
//        std::map< double, Eigen::VectorXd > interpolatedBenchmarkResult;
//        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > benchmarkInterpolator =
//                interpolators::createOneDimensionalInterpolator(
//                    propagatedStateHistory, interpolatorSettings );

//        for( auto stateIterator : propagatedHistoryBenchMark )
//        {
//            double currentTime = stateIterator.first;
//            interpolatedBenchmarkResult[ currentTime ] =
//                    benchmarkInterpolator->interpolate( currentTime );
//        }

//        input_output::writeDataMapToTextFile( interpolatedBenchmarkResult, "interpolatedStateHistory" + std::to_string(i) + "_" + behindText+ ".dat", outputPath );
//        input_output::writeDataMapToTextFile( integratorData, "integratorData" + std::to_string(i) + "_" + behindText + ".dat", outputPath );

//        input_output::writeDataMapToTextFile( propagatedStateHistory, "stateHistory" + std::to_string(i) + "_" + behindText+ ".dat", outputPath );

//        input_output::writeDataMapToTextFile( dependentVariableHistory, "dependentVariables"+ std::to_string(i) + "_" + behindText + ".dat", outputPath );
////        input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
////                                             thrustParameters ), "thrustParameters"+ std::to_string(i) + ".dat", 16, outputPath );
//    }
        // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
        return EXIT_SUCCESS;
        std::cout << "Propagation was success." << std::endl;
}

