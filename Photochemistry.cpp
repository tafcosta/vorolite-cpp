/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"


Photochemistry::Photochemistry(Mesh& mesh, double ionisationCrossSection, double recombinationCrossSection) : mesh(mesh), ionisationCrossSection(ionisationCrossSection), recombinationCrossSection(recombinationCrossSection) {
	// TODO Auto-generated constructor stub

}

void Photochemistry::evolveIonisation(double dtime) {
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

    	double dtime_in_cgs    = dtime * mesh.unitLength / mesh.unitVelocity;
        double x0              = mesh.getHIIFraction(iCell);
    	double localHIcolumn   = mesh.cellLocalColumn[iCell]  * (1 - x0) * mesh.unitMass / mesh.protonMass / mesh.unitLength / mesh.unitLength;
        double absorbedFlux    = mesh.cellIncomingFlux[iCell] * (1 - std::exp(-localHIcolumn * ionisationCrossSection));
        double nH              = mesh.getNumberDensity_in_cgs(iCell);
        double volume          = mesh.getMass(iCell)/mesh.getDensity(iCell) * mesh.unitLength * mesh.unitLength * mesh.unitLength;
        double electronDensity = mesh.getElectronNumberDensity_in_cgs(iCell);

        auto computeRate = [&](double x) {
        	return getIonisationRate(x, volume, absorbedFlux, nH) - getRecombinationRate(x, electronDensity);
        };

        double k1 = computeRate(x0);
        double k2 = computeRate(x0 + 0.5 * k1 * dtime_in_cgs);
        double k3 = computeRate(x0 + 0.5 * k2 * dtime_in_cgs);
        double k4 = computeRate(x0 + k3 * dtime_in_cgs);

        double delta = (dtime_in_cgs / 6.0) * (k1 + 2*k2 + 2*k3 + k4);

    	mesh.setHIIFraction(iCell, x0 + delta);
    }
}

double Photochemistry::getIonisationRate(double xHII, double volume, double absorbedPhotonsPerSecond, double nH){
    return absorbedPhotonsPerSecond / (nH * volume);
}


double Photochemistry::getRecombinationRate(double xHII, double electronDensity){
	return  electronDensity * xHII * recombinationCrossSection;
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}

