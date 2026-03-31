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
}

void Photochemistry::evolveIonisation(double dtime) {
    const double dtime_in_cgs = dtime * mesh.unitLength / mesh.unitVelocity;
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

    	const double oldxH = mesh.xH_old[iCell];
    	const double avgxH = mesh.xH_pred[iCell];
    	const double yHe   = mesh.getHeIIFraction(iCell);
    	const double zHe   = mesh.getHeIIIFraction(iCell);
    	const double nH    = mesh.getHNumberDensity_in_cgs(iCell);
    	const double nHe   = mesh.getHeNumberDensity_in_cgs(iCell);
    	const double temperature = mesh.getTemperature_in_K(iCell);

    	const double electronDensity = avgxH * nH + (yHe + 2.0 * zHe) * nHe;
        const double volume = mesh.getMass(iCell) / mesh.getDensity(iCell) * std::pow(mesh.unitLength, 3);

        double HydrogenNeutralFraction = std::max(1.0 - avgxH, 1e-10);
        double Gamma = 0.0;
        if (HydrogenNeutralFraction > 0 && nH > 0.0 && volume > 0.0)
            Gamma = std::max(mesh.getPhotonAbsorptionRateHI(iCell) / (HydrogenNeutralFraction * nH * volume), 0.0);

    	const double equilibriumTime = 1.0 / (Gamma + getHIIrecombinationCoefficient(temperature) * electronDensity + getHIcollisionalIonisationCoefficient(temperature) * electronDensity);
    	const double equilibriumXH   = (Gamma + getHIcollisionalIonisationCoefficient(temperature) * electronDensity) / (Gamma + getHIIrecombinationCoefficient(temperature) * electronDensity + getHIcollisionalIonisationCoefficient(temperature) * electronDensity);

    	double xHnew = equilibriumXH + (oldxH - equilibriumXH) * std::exp(-dtime_in_cgs / equilibriumTime);

    	mesh.setHIIFraction(iCell, xHnew);

    }
}


void Photochemistry::predictIonisation(double dtime)
{
    const double dtime_in_cgs = dtime * mesh.unitLength / mesh.unitVelocity;
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

    	const double xH   = mesh.xH_old[iCell];
    	const double yHe  = mesh.getHeIIFraction(iCell);
    	const double zHe  = mesh.getHeIIIFraction(iCell);
        const double nH   = mesh.getHNumberDensity_in_cgs(iCell);
        const double nHe  = mesh.getHeNumberDensity_in_cgs(iCell);
        const double temperature = mesh.getTemperature_in_K(iCell);

        const double electronDensity = xH * nH + (yHe + 2.0 * zHe) * nHe;
        const double volume = mesh.getMass(iCell) / mesh.getDensity(iCell) * std::pow(mesh.unitLength, 3);

        const double HydrogenNeutralFraction = std::max(1.0 - xH, 1e-10);
        double Gamma = 0.0;
        if (HydrogenNeutralFraction > 0 && nH > 0.0 && volume > 0.0)
        	Gamma = std::max(mesh.getPhotonAbsorptionRateHI(iCell) / (HydrogenNeutralFraction * nH * volume), 0.0);

        double equilibriumTime = 1./(Gamma + getHIIrecombinationCoefficient(temperature) * electronDensity + getHIcollisionalIonisationCoefficient(temperature) * electronDensity);
        double equilibriumXH   = (Gamma + getHIcollisionalIonisationCoefficient(temperature) * electronDensity) / (Gamma + getHIIrecombinationCoefficient(temperature) * electronDensity + getHIcollisionalIonisationCoefficient(temperature) * electronDensity);

        const double eps = dtime_in_cgs / equilibriumTime;
        const double avgFactor = (std::abs(eps) < 1e-8) ? (1.0 - 0.5*eps + eps*eps/6.0) : (-std::expm1(-eps) / eps);

        mesh.xH_pred[iCell] = equilibriumXH + (xH - equilibriumXH) * avgFactor;

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

void Photochemistry::storeOldIonisation(){
	for(int iCell = 0; iCell < mesh.numCells; ++iCell){
		mesh.xH_old[iCell] = mesh.getHIIFraction(iCell);
	}
}

double Photochemistry::solveQuadratic(double a, double b, double c, double sign)
{
// zeros of a x**2 + b x + c
	if (fabs(a) <= 1e-20){
		return -c / b;
	}
	else {
	   return solveQuadratic(b/a, c/a, 1.0 * sign);
   }
}

double Photochemistry::solveQuadratic(double p, double q, double sign){
	return -p/2 + sign * sqrt(p*p/4 - q);
}


Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}
