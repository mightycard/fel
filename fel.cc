#include <iostream>

#include "parameter.h"
#include "options.h"
#include "math.h"

int main( int argc, char ** argv )
{
    // Defining parameters and options
    Options myOptions( argc, argv, 0);

    double const gamma = myOptions.beamEnergy / parameter_e_mass;
    double const betaGamma = sqrt(gamma*gamma - 1.0);

    double sigmaX = sqrt(myOptions.bex * myOptions.emx);
    double sigmaY = sqrt(myOptions.bey * myOptions.emy);

    double const charge = myOptions.bunchCharge;
    double lineDensity = charge / gamma / myOptions.lb;

    double bc = betaGamma / gamma * parameter_c;
    double const k = 2.0 * parameter_pi * gamma / myOptions.lu;
    double const w = k * bc;

    double x0 = charge * myOptions.bMax / ( gamma * gamma * parameter_e_mass_kg
            * k * k * betaGamma * parameter_c );

    double a = 25.0 * sigmaX;
    double b = 25.0 * sigmaY;

    double const aTilda = a / parameter_pi;
    double const bTilda = b / parameter_pi;

    double tauMax = 1.0 * 2.0 * parameter_pi;
    double h = tauMax / static_cast<double >(myOptions.timeStep);

    double qmw = parameter_e / parameter_e_mass_kg / w;
    double awp = a * w / parameter_pi;
    double bwp = b * w / parameter_pi;

    double c2 = parameter_c * parameter_c;
    double wc2 = w * w / c2;

    
    // Test Particle Generations


    // Prepare transverse integrations


    // Track test particles


    // Calculate space charge fields

}
