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
	Rays(double crossSection, double maxRadius, std::vector<double> sourcePosition, double lumTotal, Mesh& mesh);
	virtual ~Rays();

	double crossSection;
	double maxRadius;
	double lumTotal;

	std::vector<double> sourcePosition;
	std::vector<double> phi, theta;

	std::vector<std::vector<double>> rayDirection;
	std::vector<std::vector<double>> rayPosition;

	std::vector<double> columnHI;
	std::vector<double> rayWeight;

	std::vector<double> distanceTravelled;

	std::vector<std::vector<double>> visitedCellColumn;
	std::vector<std::vector<int>> visitedCells;

	std::vector<bool> insideDomain;
	std::vector<bool> flagRay;

	int nRays;
	Mesh& mesh;

	void calculateRays();
	void doRadiativeTransfer(double time);
	void outputResults(std::string& ofileName);

protected:
	int startCell;
	int travelToNextCell(int iCell, int iRay, bool verbose);
	std::vector<int> rayTargetCell;

	void assignToHealpix(double L_total);
	void updateRayPosition(int iRay, double distance);
	void updateColumnAndFlux(int iRay);
	int findExitCellAndSetDistance(int iCell, int iRay, int& exitCell, double& distanceToExit, bool verbose);
	int modifyExitCellIfOnInterface(int iCell, int iRay, int exitCell, double distanceToExit, bool verbose);
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
