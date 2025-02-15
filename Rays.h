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
	Rays(int numRays, double maxRadius, std::vector<double> sourcePosition, Mesh& mesh);
	virtual ~Rays();

	int numRays;
	double maxRadius;

	std::vector<double> sourcePosition;
	std::vector<double> phi, theta;

	std::vector<std::vector<double>> rayDirection;
	std::vector<std::vector<double>> rayPosition;
	std::vector<std::vector<int>> visitedCells;


	std::vector<double> columnDensity;
	std::vector<int> numTraversedCells;
	std::vector<double> distanceTravelled;


	std::vector<bool> insideDomain;
	std::vector<bool> ignoreRay;


	Mesh& mesh;

	void doRayTracing();
	void outputResults();

protected:
	int startCell;
	int findNextCell(int iCell, int iRay, bool verbose);

private:
    void initializeDirections();
    bool warningIssued;
};

#endif /* RAYS_H_ */
