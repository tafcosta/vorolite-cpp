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
	Photochemistry(Mesh& mesh, double crossSection);
	virtual ~Photochemistry();

	Mesh& mesh;

	double crossSection;

	void evolveIonisation(double dtime);
};

#endif /* PHOTOCHEMISTRY_H_ */
