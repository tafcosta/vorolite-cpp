/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"

Photochemistry::Photochemistry(Mesh& mesh, double crossSection) : mesh(mesh), crossSection(crossSection) {
	// TODO Auto-generated constructor stub

}

void Photochemistry::evolveIonisation(double dtime) {
    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {
        double mass = mesh.cellMass[iCell];
        double flux = mesh.cellFlux[iCell];
        double x0   = mesh.cellHIIFraction[iCell];

        if (mass > 0.0) {

            auto computeAbsorbedFlux = [&](double x) {
                double columnHI = mesh.cellLocalColumn[iCell] * (1.0 - x);
                return (1.0 - std::exp(-crossSection * columnHI)) * flux;
            };

            auto computeRate = [&](double x) {
                return computeAbsorbedFlux(x) / mass;
            };

            double k1 = computeRate(x0);
            double k2 = computeRate(x0 + 0.5 * k1 * dtime);
            double k3 = computeRate(x0 + 0.5 * k2 * dtime);
            double k4 = computeRate(x0 + k3 * dtime);

            double delta = (dtime / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
            mesh.cellHIIFraction[iCell] = x0 + delta;
        }
    }
}




Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}

