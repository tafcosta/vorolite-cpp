/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"


Photochemistry::Photochemistry(Mesh& mesh, double ionisationCrossSection, double recombinationCoefficient) : mesh(mesh), HIionisationCrossSection(ionisationCrossSection), HIrecombinationCoefficient(recombinationCoefficient) {
	// TODO Auto-generated constructor stub
}

void Photochemistry::evolveIonisation(double dtime) {
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

    	if(mesh.cellLocalColumn[iCell] == 0)
    		continue;

    	double dtime_in_cgs    = dtime * mesh.unitLength / mesh.unitVelocity;

        double x0              = mesh.getHIIFraction(iCell);
        double y0              = mesh.getHeIIFraction(iCell);

    	double localColumn     = mesh.cellLocalColumn[iCell] / mesh.protonMass * (mesh.unitMass / (mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength)) * mesh.HubbleParam;
        double incomingFlux    = mesh.getIncomingFlux(iCell);
        double nH              = mesh.getHNumberDensity_in_cgs(iCell);
        double volume          = mesh.getMass(iCell)/mesh.getDensity(iCell) * (mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength) * mesh.HubbleParam * mesh.HubbleParam * mesh.HubbleParam;

        auto computeRate = [&](double x) {

        	if(x > 1)
        		x = 1;

        	double ne = x * nH;
        	double localHIcolumn = localColumn * (1 - x) * mesh.xHydrogen;
        	double flux = incomingFlux * (1 - std::exp(-localHIcolumn * HIionisationCrossSection));
        	return getIonisationRate(volume, flux, nH) - getRecombinationRate(x, ne);
        };

        double k1 = computeRate(x0);
        double k2 = computeRate(x0 + 0.5 * k1 * dtime_in_cgs);
        double k3 = computeRate(x0 + 0.5 * k2 * dtime_in_cgs);
        double k4 = computeRate(x0 + k3 * dtime_in_cgs);

        double delta = (dtime_in_cgs / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
        double new_x = x0 + delta;

        mesh.setHIIFraction(iCell, new_x);
    }
}


double Photochemistry::getIonisationRate(double volume, double flux, double nH){
    return flux / (nH * volume);
}

double Photochemistry::getRecombinationRate(double xHII, double electronDensity){
	return  electronDensity * xHII * HIrecombinationCoefficient;
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}
