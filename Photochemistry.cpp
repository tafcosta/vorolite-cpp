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

        const double NdotAbsorbedHI   = mesh.cellAbsorbedPhotonRateHI[iCell];
        const double NdotAbsorbedHeI  = mesh.cellAbsorbedPhotonRateHeI[iCell];
        const double NdotAbsorbedHeII = mesh.cellAbsorbedPhotonRateHeII[iCell];

        const double nH   = mesh.getHNumberDensity_in_cgs(iCell);
        const double nHe  = mesh.getHeNumberDensity_in_cgs(iCell);
        const double temp = mesh.getTemperature_in_K(iCell);

        const double volume =
            mesh.getMass(iCell) / mesh.getDensity(iCell) *
            (mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength) *
            mesh.HubbleParam * mesh.HubbleParam * mesh.HubbleParam;

        auto clampState = [&](double &x, double &y, double &z) {
            const double oneMinusTiny = 1.0 - 1e-20;

            if (x < 0.0) x = 0.0;
            if (x > oneMinusTiny) x = oneMinusTiny;

            if (y < 0.0) y = 0.0;
            if (z < 0.0) z = 0.0;

            const double yz = y + z;
            if (yz > oneMinusTiny) {
                const double inv = oneMinusTiny / yz;
                y *= inv;
                z *= inv;
            }
        };

        struct Rates { double dx, dy, dz; };

        auto computeRates = [&](double x, double y, double z) -> Rates {

            clampState(x, y, z);

            double ne = x * nH + (y + 2.0 * z) * nHe;
            double ionH = 0.0, ionHeI = 0.0, ionHeII = 0.0;

            if (volume > 0.0) {
                if (nH  > 0.0) ionH    = NdotAbsorbedHI   / (nH  * volume);
                if (nHe > 0.0) ionHeI  = NdotAbsorbedHeI  / (nHe * volume);
                if (nHe > 0.0) ionHeII = NdotAbsorbedHeII / (nHe * volume);
            }

            const double recH    = getRecombinationRate(Species::HI,    x, ne, temp);
            const double recHeII = getRecombinationRate(Species::HeII,  y, ne, temp);
            const double recHeIII= getRecombinationRate(Species::HeIII, z, ne, temp);

            const double C_HI    = getHIcollisionalIonisationCoefficient(temp);
            const double C_HeI   = getHeIcollisionalIonisationCoefficient(temp);
            const double C_HeII  = getHeIIcollisionalIonisationCoefficient(temp);

            const double collH   = (1.0 - x) * ne * C_HI;
            const double collHeI = (1.0 - y - z) * ne * C_HeI;
            const double collHeII= y * ne * C_HeII;

            const double dx = ionH    + collH - recH;
            const double dy = ionHeI  + collHeI  - ionHeII - recHeII + recHeIII - collHeII;
            const double dz = ionHeII + collHeII - recHeIII;

            return {dx, dy, dz};
        };

        Rates k1 = computeRates(xH, yHe, zHe);
        Rates k2 = computeRates(xH + 0.5 * k1.dx * dtime_in_cgs,
                                yHe + 0.5 * k1.dy * dtime_in_cgs,
                                zHe + 0.5 * k1.dz * dtime_in_cgs);
        Rates k3 = computeRates(xH + 0.5 * k2.dx * dtime_in_cgs,
                                yHe + 0.5 * k2.dy * dtime_in_cgs,
                                zHe + 0.5 * k2.dz * dtime_in_cgs);
        Rates k4 = computeRates(xH + k3.dx * dtime_in_cgs,
                                yHe + k3.dy * dtime_in_cgs,
                                zHe + k3.dz * dtime_in_cgs);

        xH  += (dtime_in_cgs / 6.0) * (k1.dx + 2.0*k2.dx + 2.0*k3.dx + k4.dx);
        yHe += (dtime_in_cgs / 6.0) * (k1.dy + 2.0*k2.dy + 2.0*k3.dy + k4.dy);
        zHe += (dtime_in_cgs / 6.0) * (k1.dz + 2.0*k2.dz + 2.0*k3.dz + k4.dz);

        clampState(xH, yHe, zHe);

        mesh.setHIIFraction(iCell,   xH);
        mesh.setHeIIFraction(iCell,  yHe);
        mesh.setHeIIIFraction(iCell, zHe);
    }
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
    // Hui & Gnedin (1997)
    // Returns C_HI in cm^3 s^-1

    assert(temperature > 0.0);
    temperature = std::max(temperature, 1e-20);

    const double T5 = temperature * 1e-5;
    const double sqrtT = std::sqrt(temperature);

    return 5.85e-11 * sqrtT / (1.0 + std::sqrt(T5))
         * std::exp(-157809.1 / temperature);
}

double Photochemistry::getHeIcollisionalIonisationCoefficient(double T)
{
    // Hui & Gnedin (1997)
    // Returns C_HeI in cm^3 s^-1
    assert(T > 0.0);
    T = std::max(T, 1e-20);

    const double T5 = T * 1e-5;
    const double sqrtT = std::sqrt(T);

    return 2.38e-11 * sqrtT / (1.0 + std::sqrt(T5))
         * std::exp(-285335.4 / T);
}

double Photochemistry::getHeIIcollisionalIonisationCoefficient(double T)
{
    // Hui & Gnedin (1997)
    // Returns C_HeI in cm^3 s^-1
    assert(T > 0.0);
    T = std::max(T, 1e-20);

    const double T5 = T * 1e-5;
    const double sqrtT = std::sqrt(T);

    return 5.68e-12 * sqrtT / (1.0 + std::sqrt(T5))
         * std::exp(-631515.0 / T);
}


Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}
