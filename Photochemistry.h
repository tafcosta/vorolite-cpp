/*
 * Photochemistry.h
 *
 *  Created on: 17 May 2025
 *      Author: ntc132
 */

#ifndef PHOTOCHEMISTRY_H_
#define PHOTOCHEMISTRY_H_

#include "Mesh.h"


class Photochemistry {
public:
	Photochemistry(Mesh& mesh, double recombinationCrossSection);
	virtual ~Photochemistry();

	Mesh& mesh;

	double recombinationCrossSection;

	void evolveIonisation(double dtime);
	double getIonisationRate(double xHII, double flux);
	double getRecombinationRate(double xHII, double electronDensity);
};

#endif /* PHOTOCHEMISTRY_H_ */
