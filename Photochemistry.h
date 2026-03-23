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

    double HIionisationCrossSection;
    double HeIionisationCrossSection;
    double HeIIionisationCrossSection;

    double getHIIrecombinationCoefficient(double temp);
    double getHeIIrecombinationCoefficient(double temp);
    double getHeIIIrecombinationCoefficient(double temp);
    double getHIcollisionalIonisationCoefficient(double temp);
    double getHeIcollisionalIonisationCoefficient(double temp);
    double getHeIIcollisionalIonisationCoefficient(double temp);
    double getRecombinationRate(Species species, double fraction, double electronDensity, double temperature);

    void evolveIonisation(double dtime);
    void predictIonisation(double dtime);
    void storeOldIonisation();

private:
    double solveQuadratic(double a, double b, double c, double sign);
    double solveQuadratic(double p, double q, double sign);
};

#endif /* PHOTOCHEMISTRY_H_ */
