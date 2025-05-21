/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"

Photochemistry::Photochemistry(Mesh& mesh) : mesh(mesh) {
	// TODO Auto-generated constructor stub

}

void Photochemistry::evolveIonisation(double dtime) {
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {
        double mass = mesh.cellMass[iCell];
        double transmittedFlux = mesh.cellFlux[iCell];

        if (mass > 0.0) {
            double absorbedFlux = 1.0 - transmittedFlux;
            double rate = absorbedFlux / mass;

            double x0 = mesh.cellHIIFraction[iCell];

            // RK4 steps
            double k1 = rate * dtime;
            double k2 = rate * dtime;  // Rate is constant → f(x) = constant
            double k3 = rate * dtime;
            double k4 = rate * dtime;

            // RK4 update
            double delta = (1.0 / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
            mesh.cellHIIFraction[iCell] = x0 + delta;
        }
    }
}



Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}

