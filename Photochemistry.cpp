/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"


Photochemistry::Photochemistry(Mesh& mesh, double recombinationCrossSection) : mesh(mesh), recombinationCrossSection(recombinationCrossSection) {
	// TODO Auto-generated constructor stub

}

void Photochemistry::evolveIonisation(double dtime) {
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

        double flux = mesh.cellIncomingFlux[iCell];
        double x0   = mesh.getHIIFraction(iCell);
        double electronDensity = mesh.getElectronNumberDensity_in_cgs(iCell);

        auto computeRate = [&](double x) {
        	return getIonisationRate(x, flux) - getRecombinationRate(x, electronDensity);
        };

        double k1 = computeRate(x0);
        double k2 = computeRate(x0 + 0.5 * k1 * dtime);
        double k3 = computeRate(x0 + 0.5 * k2 * dtime);
        double k4 = computeRate(x0 + k3 * dtime);

        double delta = (dtime / 6.0) * (k1 + 2*k2 + 2*k3 + k4);

    	mesh.setHIIFraction(iCell, x0 + delta);

    }
}

double Photochemistry::getIonisationRate(double xHII, double flux){
	return flux * (1.0 - xHII);
}

double Photochemistry::getRecombinationRate(double xHII, double electronDensity){
	return  electronDensity * xHII * recombinationCrossSection;
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}

