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
	Photochemistry(Mesh& mesh);
	virtual ~Photochemistry();

	Mesh& mesh;

	void evolveIonisation(double dtime);
};

#endif /* PHOTOCHEMISTRY_H_ */
