/*
 * Photochemistry.h
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#ifndef PHOTOCHEMISTRY_H_
#define PHOTOCHEMISTRY_H_

#include "Mesh.h"
#include "Rays.h"


class Photochemistry {
public:
    Photochemistry(Mesh& mesh, Rays& rays, double HIcross, double HeIcross, double HeIIcross);	virtual ~Photochemistry();

	Mesh& mesh;

    enum class Species { HI, HeII, HeIII };

    double getHIIrecombinationCoefficient(double temp);
    double getHeIIrecombinationCoefficient(double temp);
    double getHeIIIrecombinationCoefficient(double temp);
    double getHIcollisionalIonisationCoefficient(double temp);

    double HIionisationCrossSection;
    double HeIionisationCrossSection;
    double HeIIionisationCrossSection;

    void evolveIonisation(double dtime);
    void evolveIonisationWithAbsorption(double dtime, const std::vector<double>& absorbedRate);

    double getIonisationRate(double volume, double flux, double nH);
    double getRecombinationRate(Species species, double fraction, double electronDensity, double temperature);
};

#endif /* PHOTOCHEMISTRY_H_ */
