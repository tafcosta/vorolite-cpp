/*
 * Photochemistry.cpp
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#include "Photochemistry.h"


Photochemistry::Photochemistry(Mesh& mesh,
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

    	if(mesh.cellLocalColumn[iCell] == 0)
    		continue;

    	double dtime_in_cgs    = dtime * mesh.unitLength / mesh.unitVelocity;

        double xH  = mesh.getHIIFraction(iCell);
        double yHe = mesh.getHeIIFraction(iCell);
        double zHe = mesh.getHeIIIFraction(iCell);

    	double localColumn     = mesh.cellLocalColumn[iCell] / mesh.protonMass * (mesh.unitMass / (mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength)) * mesh.HubbleParam;
        double incomingFlux    = mesh.getIncomingFlux(iCell);

        double nH              = mesh.getHNumberDensity_in_cgs(iCell);
        double nHe             = mesh.getHeNumberDensity_in_cgs(iCell);

        double volume          = mesh.getMass(iCell)/mesh.getDensity(iCell) * (mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength * mesh.scaleFactor * mesh.unitLength) * mesh.HubbleParam * mesh.HubbleParam * mesh.HubbleParam;

        auto computeRateH = [&](double x) -> double {

            if (x > 1.0) x = 1.0;
            if (x < 0.0) x = 0.0;
            double ne = x * nH + (yHe + 2.0 * zHe) * nHe;
            double ion = getIonisationRate(volume, incomingFlux * sigma_HI, nH);
            double rec = getRecombinationRate(Species::HI, x, ne);
            return ion;//ion - rec;

        };

/*        auto computeRateHeI = [&](double y) -> double {
            if (y > 1.0) y = 1.0;
            if (y < 0.0) y = 0.0;
            double ne = xH * nH + (y + 2.0 * zHe) * nHe;
            double HeI_frac = 1.0 - y - zHe;
            if (HeI_frac < 0.0) HeI_frac = 0.0;
            double ion = getIonisationRate(volume, incomingFlux * sigma_HeI, nHe) * HeI_frac;
            double rec = getRecombinationRate(Species::HeII, y, ne);
            return ion - rec;
        };

        auto computeRateHeII = [&](double z) -> double {
            if (z > 1.0) z = 1.0;
            if (z < 0.0) z = 0.0;
            double ne = xH * nH + (yHe + 2.0 * z) * nHe;
            double HeII_frac = yHe;
            double ion = getIonisationRate(volume, incomingFlux * sigma_HeII, nHe) * HeII_frac;
            double rec = getRecombinationRate(Species::HeIII, z, ne);
            return ion - rec;
        };
*/

        double x1 = xH;
        //double y1 = yHe;
        //double z1 = zHe;
        double kx1 = computeRateH(x1);
        //double saved_y = yHe, saved_z = zHe;
        //yHe = y1; zHe = z1;
        //double ky1_part = computeRateHeI(y1);
        //double kz1 = computeRateHeII(z1);
        //double ky1 = ky1_part - kz1;
        //yHe = saved_y; zHe = saved_z;

        double x2 = xH + 0.5 * kx1 * dtime_in_cgs;
        //double y2 = yHe + 0.5 * ky1 * dtime_in_cgs;
        //double z2 = zHe + 0.5 * kz1 * dtime_in_cgs;
        //yHe = y2; zHe = z2;
        double kx2 = computeRateH(x2);
        //double ky2_part = computeRateHeI(y2);
        //double kz2 = computeRateHeII(z2);
        //double ky2 = ky2_part - kz2;
        //yHe = saved_y; zHe = saved_z;

        double x3 = xH + 0.5 * kx2 * dtime_in_cgs;
        //double y3 = yHe + 0.5 * ky2 * dtime_in_cgs;
        //double z3 = zHe + 0.5 * kz2 * dtime_in_cgs;
        //yHe = y3; zHe = z3;
        double kx3 = computeRateH(x3);
        //double ky3_part = computeRateHeI(y3);
        //double kz3 = computeRateHeII(z3);
        //double ky3 = ky3_part - kz3;
        //yHe = saved_y; zHe = saved_z;

        double x4 = xH + kx3 * dtime_in_cgs;
        //double y4 = yHe + ky3 * dtime_in_cgs;
        //double z4 = zHe + kz3 * dtime_in_cgs;
        //yHe = y4; zHe = z4;
        double kx4 = computeRateH(x4);
        //double ky4_part = computeRateHeI(y4);
        //double kz4 = computeRateHeII(z4);
        //double ky4 = ky4_part - kz4;
        //yHe = saved_y; zHe = saved_z;

        double delta_x = (dtime_in_cgs / 6.0) * (kx1 + 2.0*kx2 + 2.0*kx3 + kx4);
        //double delta_y = (dtime_in_cgs / 6.0) * (ky1 + 2.0*ky2 + 2.0*ky3 + ky4);
        //double delta_z = (dtime_in_cgs / 6.0) * (kz1 + 2.0*kz2 + 2.0*kz3 + kz4);

        xH  += delta_x;
        //yHe += delta_y;
        //zHe += delta_z;

        if (xH < 0.0)  xH  = 0.0; if (xH > 1.0) xH = 1.0;
        //if (yHe < 0.0) yHe = 0.0;
        //if (zHe < 0.0) zHe = 0.0;
        //double he_sum = yHe + zHe;
        //if (he_sum > 1.0) { yHe /= he_sum; zHe /= he_sum; }

        std::cout << "Volume = " << volume << ", dTime = " << dtime_in_cgs << ", scaleFactor, h = " << mesh.scaleFactor << " " << mesh.HubbleParam << std::endl;

        mesh.setHIIFraction(iCell, xH);
        mesh.setHeIIFraction(iCell, yHe);
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
