/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"


Photochemistry::Photochemistry(Mesh& mesh, double ionisationCrossSection, double recombinationCoefficient) : mesh(mesh), ionisationCrossSection(ionisationCrossSection), recombinationCoefficient(recombinationCoefficient) {
	// TODO Auto-generated constructor stub

}

void Photochemistry::evolveIonisation(double dtime) {
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

    	double dtime_in_cgs    = dtime * mesh.unitLength / mesh.unitVelocity;
        double x0              = mesh.getHIIFraction(iCell);
    	double localColumn     = mesh.cellLocalColumn[iCell] * mesh.unitMass / mesh.protonMass / mesh.unitLength / mesh.unitLength;
        double flux            = mesh.cellIncomingFlux[iCell];
        double nH              = mesh.getNumberDensity_in_cgs(iCell);
        double volume          = mesh.getMass(iCell)/mesh.getDensity(iCell) * mesh.unitLength * mesh.unitLength * mesh.unitLength;

        auto computeRate = [&](double x) {
        	double ne = x * nH;
        	return getIonisationRate(x, volume, flux, nH, localColumn) - getRecombinationRate(x, ne);
        };

        double k1 = computeRate(x0);
        double k2 = computeRate(x0 + 0.5 * k1 * dtime_in_cgs);
        double k3 = computeRate(x0 + 0.5 * k2 * dtime_in_cgs);
        double k4 = computeRate(x0 + k3 * dtime_in_cgs);

        double delta = (dtime_in_cgs / 6.0) * (k1 + 2*k2 + 2*k3 + k4);

    	mesh.setHIIFraction(iCell, x0 + delta);

    	//if(iCell == 915)
    	//	std::cout << "xHII = " << mesh.getHIIFraction(iCell) << std::endl;

    }
}

double Photochemistry::getIonisationRate(double xHII, double volume, double flux, double nH, double localColumn){

	if(xHII > 1)
		xHII = 1;

	double localHIcolumn = localColumn * (1 - xHII);
    return flux * (1 - std::exp(-localHIcolumn * ionisationCrossSection)) / (nH * volume);
}

double Photochemistry::getRecombinationRate(double xHII, double electronDensity){
	return  electronDensity * xHII * recombinationCoefficient;
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}


/*
Photochemistry::Photochemistry(Mesh& mesh, double ionisationCrossSection, double recombinationCrossSection)
    : mesh(mesh), ionisationCrossSection(ionisationCrossSection), recombinationCrossSection(recombinationCrossSection) {}

Photochemistry::~Photochemistry() {}

void Photochemistry::evolveIonisation(double dtime) {
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

        double dtime_in_cgs = dtime * mesh.unitLength / mesh.unitVelocity;
        double x0 = mesh.getHIIFraction(iCell);

        double localHIcolumn = mesh.cellLocalColumn[iCell] * (1 - x0) * mesh.unitMass / mesh.protonMass
                             / (mesh.unitLength * mesh.unitLength);
        double absorbedFlux = mesh.cellIncomingFlux[iCell] * (1 - std::exp(-localHIcolumn * ionisationCrossSection));
        double nH = mesh.getNumberDensity_in_cgs(iCell);
        double volume = mesh.getMass(iCell) / mesh.getDensity(iCell) * std::pow(mesh.unitLength, 3);
        double electronDensity = mesh.getElectronNumberDensity_in_cgs(iCell);

        auto ionisationRate = [&](double x) {
            return getIonisationRate(x, volume, absorbedFlux, nH);
        };

        auto recombinationRate = [&](double x) {
            return getRecombinationRate(x, electronDensity);
        };

        auto residual = [&](double x_next) {
            double x_mid = 0.5 * (x0 + x_next);
            return x_next - x0 - dtime_in_cgs * (ionisationRate(x_mid) - recombinationRate(x_mid));
        };

        auto residual_derivative = [&](double x_next) {
            double eps = 1e-6;
            return (residual(x_next + eps) - residual(x_next - eps)) / (2 * eps);
        };

        // Newton-Raphson iteration
        double x_next = x0;
        const int max_iter = 20;
        const double tol = 1e-8;
        for (int iter = 0; iter < max_iter; ++iter) {
            double r = residual(x_next);
            double dr = residual_derivative(x_next);
            if (std::abs(dr) < 1e-10) break;  // Avoid division by near-zero
            double dx = -r / dr;
            x_next += dx;
            if (std::abs(dx) < tol) break;
        }

        mesh.setHIIFraction(iCell, x_next);
    }
}

double Photochemistry::getIonisationRate(double xHII, double volume, double absorbedPhotonsPerSecond, double nH) {
    return absorbedPhotonsPerSecond / (nH * volume);
}

double Photochemistry::getRecombinationRate(double xHII, double electronDensity) {
    return electronDensity * xHII * recombinationCrossSection;
}
*/


