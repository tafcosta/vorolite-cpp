/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: Tiago Costa
 */

#include "Photochemistry.h"


Photochemistry::Photochemistry(Mesh& mesh, Rays& rays,
                               double HIcross, double HeIcross, double HeIIcross)
	: mesh(mesh),
      HIionisationCrossSection(HIcross),
      HeIionisationCrossSection(HeIcross),
      HeIIionisationCrossSection(HeIIcross)
{
	// TODO Auto-generated constructor stub
}

void Photochemistry::evolveIonisation(double dtime) {

    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

        const double dtime_in_cgs = dtime * mesh.unitLength / mesh.unitVelocity;

        double xH  = mesh.getHIIFraction(iCell);
        double yHe = mesh.getHeIIFraction(iCell);
        double zHe = mesh.getHeIIIFraction(iCell);

        const double NdotAbsorbed = mesh.cellAbsorbedPhotonRate[iCell];

        const double nH   = mesh.getHNumberDensity_in_cgs(iCell);
        const double nHe  = mesh.getHeNumberDensity_in_cgs(iCell);
        const double temp = mesh.getTemperature_in_K(iCell);

        const double volume =
            mesh.getMass(iCell) / mesh.getDensity(iCell) *
            (mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength) * mesh.HubbleParam * mesh.HubbleParam * mesh.HubbleParam;


        auto computeRateH = [&](double x) -> double {

            double ion = 0.0;

            const double xMax = 1.0 - 1.e-20;
            if (x > xMax) x = xMax;
            if (x < 0.0)  x = 0.0;

            double ne = x * nH + (yHe + 2.0 * zHe) * nHe;

            if (nH > 0.0 && volume > 0.0)
                ion = NdotAbsorbed / (nH * volume);  // [1/s]

            const double rec = getRecombinationRate(Species::HI, x, ne, temp);

            const double C_HI = getHIcollisionalIonisationCoefficient(temp);
            const double coll = (1.0 - x) * ne * C_HI;

            return ion + coll - rec;
        };

        double kx1 = computeRateH(xH);
        double x2  = xH + 0.5 * kx1 * dtime_in_cgs;
        double kx2 = computeRateH(x2);
        double x3  = xH + 0.5 * kx2 * dtime_in_cgs;
        double kx3 = computeRateH(x3);
        double x4  = xH + kx3 * dtime_in_cgs;
        double kx4 = computeRateH(x4);
        double delta_x = (dtime_in_cgs / 6.0) * (kx1 + 2.0*kx2 + 2.0*kx3 + kx4);

        mesh.cellNetIonisationRate[iCell] = delta_x * nH * volume / dtime_in_cgs;

        xH  += delta_x;

        if (xH < 0.0)  xH = 0.0;
    	const double xMax = 1.0 - 1.e-20;
    	if (xH > xMax) xH = xMax;

        mesh.setHIIFraction(iCell,   xH);
        mesh.setHeIIFraction(iCell,  yHe);
        mesh.setHeIIIFraction(iCell, zHe);
    }

}


double Photochemistry::getIonisationRate(double volume, double flux, double nH){
    return flux / (nH * volume);
}


double Photochemistry::getRecombinationRate(Species species, double fraction, double electronDensity, double temp) {
    double alpha = 0.0;

    switch(species) {
        case Species::HI:
            alpha = getHIIrecombinationCoefficient(temp);
            break;
        case Species::HeII:
            alpha = getHeIIrecombinationCoefficient(temp);
            break;
        case Species::HeIII:
            alpha = getHeIIIrecombinationCoefficient(temp);
            break;
    }

    return fraction * electronDensity * alpha;
}

double Photochemistry::getHIIrecombinationCoefficient(double temperature)
{
    // Hui & Gnedin (1997)
    // Returns alpha_B in cm^3 s^-1
    assert(temperature > 0.0);

    temperature = std::max(temperature, 1e-20);

    const double lambda = 315614.0 / temperature;

    return 2.753e-14
         * lambda * std::sqrt(lambda)
         * std::pow(1.0 + std::pow(lambda / 2.740, 0.407), -2.242);
}

double Photochemistry::getHeIIrecombinationCoefficient(double temperature)
{
    // Hui & Gnedin (1997)
    // Returns alpha_B in cm^3 s^-1

    assert(temperature > 0.0);
    temperature = std::max(temperature, 1e-20);

    const double lambda_HeI = 570670.0 / temperature;
    return 1.26e-14 * std::pow(lambda_HeI, 0.750);
}

double Photochemistry::getHeIIIrecombinationCoefficient(double temperature)
{
    // Hui & Gnedin (1997)
    // Returns alpha_B in cm^3 s^-1

    assert(temperature > 0.0);
    temperature = std::max(temperature, 1e-20);

    const double lambda_HeII = 1263030.0 / temperature;

    return (2.0 * 2.753e-14)
         * lambda_HeII * std::sqrt(lambda_HeII)
         * std::pow(1.0 + std::pow(lambda_HeII / 2.740, 0.407), -2.242);
}


double Photochemistry::getHIcollisionalIonisationCoefficient(double temperature)
{
    // Returns C_HI(temperature) in cm^3 s^-1

    assert(temperature > 0.0);
    temperature = std::max(temperature, 1e-20);

    const double T5 = temperature * 1e-5;
    const double sqrtT = std::sqrt(temperature);

    return 5.85e-11 * sqrtT / (1.0 + std::sqrt(T5))
         * std::exp(-157809.1 / temperature);
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}
