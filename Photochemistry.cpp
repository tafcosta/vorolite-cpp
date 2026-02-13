/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"


Photochemistry::Photochemistry(Mesh& mesh, Rays& rays,
                               double HIcross, double HIrecomb,
                               double HeIcross, double HeIrecomb,
                               double HeIIcross, double HeIIrecomb)
	: mesh(mesh),
      HIionisationCrossSection(HIcross),
      HIrecombinationCoefficient(HIrecomb),
      HeIionisationCrossSection(HeIcross),
      HeIrecombinationCoefficient(HeIrecomb),
      HeIIionisationCrossSection(HeIIcross),
      HeIIrecombinationCoefficient(HeIIrecomb)
{
	// TODO Auto-generated constructor stub
}

void Photochemistry::evolveIonisation(double dtime) {

    const double sigma_HI   = HIionisationCrossSection;
    const double sigma_HeI  = HeIionisationCrossSection;
    const double sigma_HeII = HeIIionisationCrossSection;

    for (int iCell = 0; iCell < mesh.numCells; ++iCell) {

        const double dtime_in_cgs = dtime * mesh.unitLength / mesh.unitVelocity;

        double xH  = mesh.getHIIFraction(iCell);
        double yHe = mesh.getHeIIFraction(iCell);
        double zHe = mesh.getHeIIIFraction(iCell);

        const double NdotAbsorbed = mesh.cellAbsorbedPhotonRate[iCell];

        const double nH  = mesh.getHNumberDensity_in_cgs(iCell);   // [cm^-3]
        const double nHe = mesh.getHeNumberDensity_in_cgs(iCell);  // [cm^-3]

        const double volume =
            mesh.getMass(iCell) / mesh.getDensity(iCell) *
            (mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength) * mesh.HubbleParam * mesh.HubbleParam * mesh.HubbleParam;

        if(iCell == 8820){
        	std::cout << "nAbsCand = " << NdotAbsorbed * dtime_in_cgs << ", nNeutrals = " << nH * volume * (1- xH) << std::endl;
        }


        auto computeRateH = [&](double x) -> double {

        	double ion = 0.;
        	const double xMax = 1.0 - 1.e-20;
        	if (x > xMax) x = xMax;
            if (x < 0.0)  x = 0.0;

            double ne = x * nH + (yHe + 2.0 * zHe) * nHe;
            double ion_rate = 0.0;

            if (nH > 0.0 && volume > 0.0)
                ion = NdotAbsorbed / (nH * volume);  // [1/s]


            double rec = getRecombinationRate(Species::HI, x, ne);
            return ion-rec;

        };

        double kx1 = computeRateH(xH);
        double x2  = xH + 0.5 * kx1 * dtime_in_cgs;
        double kx2 = computeRateH(x2);
        double x3  = xH + 0.5 * kx2 * dtime_in_cgs;
        double kx3 = computeRateH(x3);
        double x4  = xH + kx3 * dtime_in_cgs;
        double kx4 = computeRateH(x4);
        double delta_x = (dtime_in_cgs / 6.0) * (kx1 + 2.0*kx2 + 2.0*kx3 + kx4);

        mesh.cellNetIonisationRate[iCell] = delta_x * nH * volume / dtime_in_cgs;

        xH  += delta_x;

        if (xH < 0.0)  xH = 0.0;
    	const double xMax = 1.0 - 1.e-20;
    	if (xH > xMax) xH = xMax;

        mesh.setHIIFraction(iCell,   xH);
        mesh.setHeIIFraction(iCell,  yHe);
        mesh.setHeIIIFraction(iCell, zHe);
    }

}


double Photochemistry::getIonisationRate(double volume, double flux, double nH){
    return flux / (nH * volume);
}


double Photochemistry::getRecombinationRate(Species species, double fraction, double electronDensity) {
    double alpha = 0.0;

    switch(species) {
        case Species::HI:
            alpha = HIrecombinationCoefficient;
            break;
        case Species::HeII:
            alpha = HeIrecombinationCoefficient;
            break;
        case Species::HeIII:
            alpha = HeIIrecombinationCoefficient;
            break;
    }

    return fraction * electronDensity * alpha;
}

Photochemistry::~Photochemistry() {
	// TODO Auto-generated destructor stub
}
