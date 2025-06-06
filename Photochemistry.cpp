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

    	if(mesh.cellLocalColumn[iCell] == 0)
    		continue;

    	double dtime_in_cgs    = dtime * mesh.unitLength / mesh.unitVelocity;
        double x0              = mesh.getHIIFraction(iCell);
    	double localColumn     = mesh.cellLocalColumn[iCell] * mesh.unitMass / mesh.protonMass / mesh.unitLength / mesh.unitLength;
        double incomingFlux    = mesh.cellIncomingFlux[iCell];//mesh.getFlux(iCell);
        double nH              = mesh.getNumberDensity_in_cgs(iCell);
        double volume          = mesh.getMass(iCell)/mesh.getDensity(iCell) * mesh.unitLength * mesh.unitLength * mesh.unitLength;

        std::cout << incomingFlux << std::endl;

        auto computeRate = [&](double x) {

        	if(x > 1)
        		x = 1;

        	double ne = x * nH;
        	double localHIcolumn = localColumn * (1 - x);
        	double flux = incomingFlux * (1 - std::exp(-localHIcolumn * ionisationCrossSection));
        	return getIonisationRate(volume, flux, nH) - getRecombinationRate(x, ne);
        };

        double k1 = computeRate(x0);
        double k2 = computeRate(x0 + 0.5 * k1 * dtime_in_cgs);
        double k3 = computeRate(x0 + 0.5 * k2 * dtime_in_cgs);
        double k4 = computeRate(x0 + k3 * dtime_in_cgs);

        double delta = (dtime_in_cgs / 6.0) * (k1 + 2*k2 + 2*k3 + k4);

        /*
        if(iCell == 432){
        	double xnew = x0 + delta;
        	if(xnew > 1)
        		xnew = 1;

        	double localHIcolumn = localColumn * (1 - xnew);
        	std::cout << xnew << " " << getIonisationRate(volume, incomingFlux * (1 - std::exp(-localHIcolumn * ionisationCrossSection)), nH) << " " << getRecombinationRate(xnew, xnew * nH) << std::endl;
        }
        */

    	mesh.setHIIFraction(iCell, x0 + delta);

    }
}


double Photochemistry::getIonisationRate(double volume, double flux, double nH){
    return flux / (nH * volume);
}

double Photochemistry::getRecombinationRate(double xHII, double electronDensity){
	return  electronDensity * xHII * recombinationCoefficient;
}


Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}
