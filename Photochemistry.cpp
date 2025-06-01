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

        double x0 = mesh.getHIIFraction(iCell);
    	double localColumn = mesh.cellLocalColumn[iCell] * mesh.unitMass/mesh.unitLength/mesh.unitLength / mesh.protonMass;
        double absorbedFlux = mesh.cellIncomingFlux[iCell]  * (1 - std::exp(-localColumn * ionisationCrossSection));
        double nH = mesh.getNumberDensity_in_cgs(iCell);
        double volume = mesh.getMass(iCell)/mesh.getDensity(iCell);
        double electronDensity = mesh.getElectronNumberDensity_in_cgs(iCell);

        std::cout << iCell << " " << mesh.cellLocalColumn[iCell] << std::endl;


        auto computeRate = [&](double x) {
        	return getIonisationRate(x, volume, absorbedFlux, nH) - getRecombinationRate(x, electronDensity);
        };

        double k1 = computeRate(x0);
        double k2 = computeRate(x0 + 0.5 * k1 * dtime);
        double k3 = computeRate(x0 + 0.5 * k2 * dtime);
        double k4 = computeRate(x0 + k3 * dtime);

        double delta = (dtime / 6.0) * (k1 + 2*k2 + 2*k3 + k4);

    	mesh.setHIIFraction(iCell, x0 + delta);

    }
}

double Photochemistry::getIonisationRate(double xHII, double volume, double absorbedPhotonsPerSecond, double nH){
    double nHI = (1.0 - xHII) * nH;
    double N_HI = nHI * volume;


    return absorbedPhotonsPerSecond / N_HI; // per-atom ionization rate
}


double Photochemistry::getRecombinationRate(double xHII, double electronDensity){
	return  electronDensity * xHII * recombinationCrossSection;
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}

