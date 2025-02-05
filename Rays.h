/*
 * Rays.h
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#ifndef RAYS_H_
#define RAYS_H_

#include "Mesh.h"

class Rays {
public:
	Rays(int numRays, std::vector<double> sourcePosition, Mesh& mesh);
	virtual ~Rays();

	int numRays;
	std::vector<double> sourcePosition;
	std::vector<double> phi, theta;

	std::vector<std::vector<double>> rayDirection;
	std::vector<std::vector<double>> rayPosition;

	std::vector<double> opticalDepth;
	std::vector<double> numTraversedCells;
	std::vector<bool> insideDomain;

	Mesh& mesh;

	void doRayTracing();

protected:
	int startCell;
	int findNextCell(int iCell, int iRay);

private:
    void initializeDirections();
};

#endif /* RAYS_H_ */
