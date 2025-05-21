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
	Rays(double maxRadius, std::vector<double> sourcePosition, Mesh& mesh);
	virtual ~Rays();

	double maxRadius;

	std::vector<double> sourcePosition;
	std::vector<double> phi, theta;

	std::vector<std::vector<double>> rayDirection;
	std::vector<std::vector<double>> rayPosition;

	std::vector<double> columnHI;

	std::vector<int> lastVisitedCell;
	std::vector<int> numTraversedCells;
	std::vector<double> distanceTravelled;

	std::vector<bool> insideDomain;
	std::vector<bool> flagRay;

	int nRays;
	Mesh& mesh;

	void doRayTracing();
	void outputResults(std::string& ofileName);

protected:
	int startCell;
	int travelToNextCell(int iCell, int iRay, bool verbose);
	std::vector<int> rayFinalCell;

	void findExitCellAndSetDistance(int iCell, int iRay, int& exitCell, double& exitToCell, bool verbose);
	int modifyExitCellIfOnInterface(int iCell, int iRay, int exitCell, double distanceToExit, bool verbose);
	bool shouldRayBeTerminated(int iRay, double distanceToExit);
	bool updateRayAndIsMaxReached(int iCell, int iRay, double& distanceToExit);
	double getOvershootDistance(int exitCell, int iRay, double distanceToExit, bool verbose);

private:
    void initializeDirections();
    void initializePositions();

    void setNumRays();
    bool warningIssued;

    double minTolerance = 1.e-7;
    int maxnIter = 100;
};

#endif /* RAYS_H_ */
