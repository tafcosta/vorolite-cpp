/*
 * Rays.h
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#ifndef RAYS_H_
#define RAYS_H_

#include "Mesh.h"
#include "Source.h"

class Rays {
public:
	Rays(double ionisationCrossSectionHI, double ionisationCrossSectionHeI, double ionisationCrossSectionHeII, double maxRadius, int64_t Nside, Mesh& mesh, Source& source);
	virtual ~Rays();

	bool timeDependent = false;

	double ionisationCrossSectionHI;
	double ionisationCrossSectionHeI;
	double ionisationCrossSectionHeII;

	double ionisationCrossSectionHI_inInternalUnits;
	double ionisationCrossSectionHeI_inInternalUnits;
	double ionisationCrossSectionHeII_inInternalUnits;

	double dustAbsorptionOpacity;
	double dustAbsorptionOpacity_inInternalUnits;
	double maxRadius;
	int64_t Nside;

	std::vector<double> sourcePosition;
	std::vector<double> phi, theta;

	std::vector<std::vector<double>> rayDirection;
	std::vector<std::vector<double>> rayPosition;

	std::vector<double> columnHI;
	std::vector<double> columnDust;

	std::vector<double> rayWeight;
	std::vector<double> distanceTravelled;
	std::vector<double> finalLuminosity;


	std::vector<std::vector<double>> visitedCellColumn;
	std::vector<std::vector<int>> visitedCells;

	std::vector<bool> insideDomain;
	std::vector<bool> flagRay;

	int nRays;
	Mesh& mesh;
	Source& source;

	const double speedOfLight = 2.99792458e10;
	double speedOfLightInternal = speedOfLight/mesh.unitVelocity;

	void calculateRays();
	void doRadiativeTransfer(double time, double dtime);
	void outputResults(std::string& ofileName);

protected:
	int startCell;
	int travelToNextCell(int iCell, int iRay, bool verbose);
	std::vector<int> rayTargetCell;

	void assignToHealpix(int64_t healpixNside);
	void updateRayPosition(int iRay, double distance);
	void updateColumnAndFlux(int iRay, double time, double dtime);
	int findExitCellAndSetDistance(int iCell, int iRay, int& exitCell, double& distanceToExit, bool verbose);
	int modifyExitCellIfOnInterface(int iCell, int iRay, int exitCell, double distanceToExit, bool verbose);
	bool updateRayAndIsMaxReached(int iCell, int iRay, double& distanceToExit);
	double getOvershootDistance(int exitCell, int iRay, double distanceToExit, bool verbose);
	double distanceSquared(std::vector<float>& a, std::vector<float>& b);

private:
    void initializeDirections();
    void initializePositions();

    void setNumRays();
    bool warningIssued;

    double minTolerance = 1.e-7;
    int maxnIter = 100;
};

#endif /* RAYS_H_ */
