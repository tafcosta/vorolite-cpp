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

    	const double x0   = mesh.xH_old[iCell];
    	const double xAvg = mesh.xH_pred[iCell];

    	const double yHe  = mesh.getHeIIFraction(iCell);
    	const double zHe  = mesh.getHeIIIFraction(iCell);

    	const double nH   = mesh.getHNumberDensity_in_cgs(iCell);
    	const double nHe  = mesh.getHeNumberDensity_in_cgs(iCell);

    	const double temp = mesh.getTemperature_in_K(iCell);
    	const double ne   = xAvg * nH + (yHe + 2.0 * zHe) * nHe;

        const double volume =
            mesh.getMass(iCell) / mesh.getDensity(iCell) *
            std::pow(mesh.scaleFactor * mesh.unitLength, 3) *
            std::pow(mesh.HubbleParam, 3);

        double neutral = std::max(1.0 - xAvg, 1e-8);
        double Gamma = 1.e-20;
        if (nH > 0.0 && volume > 0.0)
            Gamma = std::max(mesh.cellPhotonAbsorptionRateHI[iCell] / (neutral * nH * volume), 1.e-20);

    	const double alpha = getHIIrecombinationCoefficient(temp);
    	const double equilibriumTime = 1.0 / (Gamma + alpha * ne);
    	const double equilibriumXH   = Gamma / (Gamma + alpha * ne);

    	double xHnew = equilibriumXH + (x0 - equilibriumXH) * std::exp(-dtime_in_cgs / equilibriumTime);

    	mesh.setHIIFraction(iCell, xHnew);

    	/*
    	if(iCell == 16685){

        	std::cout << nH << " " << mesh.xH_old[iCell] << " " << mesh.xH_pred[iCell] << " " << xHnew << " " << Gamma << std::endl;

        }*/

        mesh.setHIIFraction(iCell,  xHnew);

    }
}


/*
void Photochemistry::evolveIonisation(double dtime) {

    const double dtime_in_cgs = dtime * mesh.unitLength / mesh.unitVelocity;

    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {
        double xH  = mesh.xH_pred[iCell];
        double yHe = mesh.getHeIIFraction(iCell);
        double zHe = mesh.getHeIIIFraction(iCell);

        const double flux_in_cgs = mesh.cellFlux[iCell];
        const double nH     = mesh.getHNumberDensity_in_cgs(iCell);
        const double nHe    = mesh.getHeNumberDensity_in_cgs(iCell);
        const double temp   = mesh.getTemperature_in_K(iCell);
        const double volume = mesh.getMass(iCell) / mesh.getDensity(iCell) *
            std::pow(mesh.scaleFactor * mesh.unitLength, 3) * std::pow(mesh.HubbleParam, 3);

        struct Rates { double dx, dy, dz; };

        auto computeRates = [&](double x, double y, double z) -> Rates {
            const double ne = x * nH + (y + 2.0 * z) * nHe;

            double ionH = 0.0, ionHeI = 0.0, ionHeII = 0.0;
            if (volume > 0.0) {

            	if (nH > 0.0) {
            		ionH = mesh.cellPhotonAbsorptionRateHI[iCell] / (nH * volume);
            	}

                if (nHe > 0.0) ionHeI  = (1.0 - y) * flux_in_cgs * 0.;
                if (nHe > 0.0) ionHeII = (1.0 - z) * flux_in_cgs * 0.;
            }

            const double recH     = getRecombinationRate(Species::HI,    x, ne, temp);
            const double recHeII  = getRecombinationRate(Species::HeII,  y, ne, temp);
            const double recHeIII = getRecombinationRate(Species::HeIII, z, ne, temp);

            const double C_HI     = 0.0;//getHIcollisionalIonisationCoefficient(temp);
            const double C_HeI    = 0.0;//getHeIcollisionalIonisationCoefficient(temp);
            const double C_HeII   = 0.0;//getHeIIcollisionalIonisationCoefficient(temp);

            const double collH    = (1.0 - x) * ne * C_HI;
            const double collHeI  = (1.0 - y - z) * ne * C_HeI;
            const double collHeII = y * ne * C_HeII;

            const double dx = ionH    + collH    - recH;
            const double dy = ionHeI  + collHeI  - ionHeII - recHeII + recHeIII - collHeII;
            const double dz = ionHeII + collHeII - recHeIII;

            return {dx, dy, dz};
        };

        Rates k1 = computeRates(xH, yHe, zHe);

        Rates k2 = computeRates(xH  + 0.5 * k1.dx * dtime_in_cgs,
        		yHe + 0.5 * k1.dy * dtime_in_cgs,
				zHe + 0.5 * k1.dz * dtime_in_cgs);

        Rates k3 = computeRates(xH  + 0.5 * k2.dx * dtime_in_cgs,
        		yHe + 0.5 * k2.dy * dtime_in_cgs,
				zHe + 0.5 * k2.dz * dtime_in_cgs);

        Rates k4 = computeRates(xH + k3.dx * dtime_in_cgs,
        		yHe + k3.dy * dtime_in_cgs,
				zHe + k3.dz * dtime_in_cgs);


        xH  += (dtime_in_cgs / 6.0) * (k1.dx + 2.0*k2.dx + 2.0*k3.dx + k4.dx);
        yHe += (dtime_in_cgs / 6.0) * (k1.dy + 2.0*k2.dy + 2.0*k3.dy + k4.dy);
        zHe += (dtime_in_cgs / 6.0) * (k1.dz + 2.0*k2.dz + 2.0*k3.dz + k4.dz);

        xH = std::clamp(xH, 0.0, 1.0);

        mesh.setHIIFraction(iCell,   xH);
        mesh.setHeIIFraction(iCell,  yHe);
        mesh.setHeIIIFraction(iCell, zHe);
    }
}*/


void Photochemistry::predictIonisation(double dtime)
{
    const double dtime_in_cgs = dtime * mesh.unitLength / mesh.unitVelocity;
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

    	const double xH   = mesh.xH_old[iCell];
    	const double yHe  = mesh.getHeIIFraction(iCell);
    	const double zHe  = mesh.getHeIIIFraction(iCell);

        const double nH   = mesh.getHNumberDensity_in_cgs(iCell);
        const double nHe  = mesh.getHeNumberDensity_in_cgs(iCell);

        const double temp = mesh.getTemperature_in_K(iCell);
        const double ne   = xH * nH + (yHe + 2.0 * zHe) * nHe;

        const double volume =
            mesh.getMass(iCell) / mesh.getDensity(iCell) *
            std::pow(mesh.scaleFactor * mesh.unitLength, 3) *
            std::pow(mesh.HubbleParam, 3);

        double neutral = std::max(1.0 - xH, 1e-8);
        double Gamma = 1.e-20;
        if (nH > 0.0 && volume > 0.0)
            Gamma = std::max(mesh.cellPhotonAbsorptionRateHI[iCell] / (neutral * nH * volume), 1.e-20);

        double equilibriumTime = 1./(Gamma + getHIIrecombinationCoefficient(temp) * ne);
        double equilibriumXH   = Gamma / (Gamma + getHIIrecombinationCoefficient(temp) * ne);

        mesh.xH_pred[iCell] = equilibriumXH + (xH - equilibriumXH) * (1 - std::exp(-dtime_in_cgs/equilibriumTime)) * equilibriumTime/dtime_in_cgs;

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

void Photochemistry::storeOldIonisation(){
	for(int iCell = 0; iCell < mesh.numCells; ++iCell){
		mesh.xH_old[iCell] = mesh.getHIIFraction(iCell);
	}
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}
