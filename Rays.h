/*
 * Rays.h
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#ifndef RAYS_H_
#define RAYS_H_

class Rays {
public:
	Rays(int numRays, std::vector<double> sourcePosition);
	virtual ~Rays();

	int numRays;
	std::vector<double> sourcePosition;
	std::vector<double> phi, theta;

	std::vector<std::vector<double>> direction;
	std::vector<std::vector<double>> position;

	std::vector<double> opticalDepth;
	std::vector<double> numTraversedCells;
	std::vector<bool> insideDomain;

	void doRayTracing();

private:
    void initializeDirections();
};

#endif /* RAYS_H_ */
