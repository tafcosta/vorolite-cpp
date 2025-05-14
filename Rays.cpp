/*
 * Rays.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Rays.h"
#include "Mesh.h"

Rays::Rays(double maxRadius, std::vector<double> sourcePosition, Mesh& mesh) : maxRadius(maxRadius), sourcePosition(sourcePosition), mesh(mesh) {
	startCell = mesh.findHostCellID(sourcePosition, -1)[0];

	setNumRays();

	rayFinalCell = std::vector<int> (nRays);
    for (int iRay = 0; iRay < nRays; ++iRay)
    	rayFinalCell[iRay] = iRay;

	rayDirection = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	theta = std::vector<double>(nRays, 0.0);
	phi   = std::vector<double>(nRays, 0.0);
	initializeDirections();

	rayPosition = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	initializePositions();

	columnDensity     = std::vector<double>(nRays, 0.0);

	distanceTravelled = std::vector<double>(nRays, 0.0);
	insideDomain      = std::vector<bool>(nRays, true);
	flagRay           = std::vector<bool>(nRays, false);
	numTraversedCells = std::vector<int>(nRays, 0.0);
	lastVisitedCell   = std::vector<int>(nRays, 0.0);

	warningIssued = false;
}

void Rays::setNumRays(){

	nRays = mesh.numCells;
}

void Rays::initializeDirections() {
	std::vector<float> cellPos;
    double xDistance, yDistance, zDistance;
    double rDistance;

    for (int iRay = 0; iRay < nRays; ++iRay) {

    	cellPos   = mesh.cellCoordinates[iRay];
    	xDistance = cellPos[0] - mesh.cellCoordinates[startCell][0];//sourcePosition[0];
    	yDistance = cellPos[1] - mesh.cellCoordinates[startCell][1];//sourcePosition[1];
    	zDistance = cellPos[2] - mesh.cellCoordinates[startCell][2];//sourcePosition[2];
    	rDistance = std::sqrt(xDistance*xDistance + yDistance*yDistance + zDistance*zDistance);

        rayDirection[iRay][0] = xDistance / rDistance;
        rayDirection[iRay][1] = yDistance / rDistance;
        rayDirection[iRay][2] = zDistance / rDistance;

        if(iRay == 9261){
        	std::cout << "DEBUG: Setting Ray Dir = " << rayDirection[iRay][0] << " " << rayDirection[iRay][1] << " " << rayDirection[iRay][2] << std::endl;
        }

        phi[iRay]   = std::atan2(yDistance, xDistance);
        theta[iRay] = std::acos(zDistance / rDistance);
    }
}

void Rays::initializePositions() {
    for (int iRay = 0; iRay < nRays; ++iRay)
        for (int i = 0; i < 3; ++i)
        	rayPosition[iRay][i] = 	mesh.cellCoordinates[startCell][i];
}


 int Rays::travelToNextCell(int iCell, int iRay, bool verbose){
	 double distanceToExit = std::numeric_limits<double>::max();
	 double overshoot = 0.;
	 int exitCell = -1;

	 if(verbose){
		 std::vector<float> cellPos = mesh.cellCoordinates[iCell];
		 double distanceBetRayAndCell = sqrt((cellPos[0] - rayPosition[iRay][0]) * (cellPos[0] - rayPosition[iRay][0]) + (cellPos[1] - rayPosition[iRay][1]) * (cellPos[1] - rayPosition[iRay][1]) + (cellPos[2] - rayPosition[iRay][2]) * (cellPos[2] - rayPosition[iRay][2]));
	 
		 std::cout << "Host Index = " << iCell << " Host Cell Position = " << cellPos[0] << ", " << cellPos[1] << ", " << cellPos[2] << " Distance from Ray to Cell (before update) = " << distanceBetRayAndCell << "\n";
	 }


	 findExitCellAndSetDistance(iCell, iRay, exitCell, distanceToExit, verbose);
	 exitCell = modifyExitCellIfOnInterface(iCell, iRay, exitCell, distanceToExit, verbose);

	 if(shouldRayBeTerminated(iRay, distanceToExit))
		 insideDomain[iRay] = false;

	 if(insideDomain[iRay]){

		 if(updateRayAndIsMaxReached(iCell, iRay, distanceToExit))
			 insideDomain[iRay] = false;

		 overshoot = getOvershootDistance(exitCell, iRay, distanceToExit, verbose);

		 if(updateRayAndIsMaxReached(exitCell, iRay, overshoot))
			 insideDomain[iRay] = false;
	 }

	 for (int i = 0; i < 3; i++)
		 rayPosition[iRay][i] += rayDirection[iRay][i] * (distanceToExit + overshoot);

	 if(insideDomain[iRay] == false)
		 return -1;





    if(verbose){
    	double distanceToCell = mesh.getDistanceToCell(rayPosition[iRay], exitCell);
    	std::vector<int> closestCells;

    	std::cout << "Ray position (after update) = " << rayPosition[iRay][0] << " "
    			<< rayPosition[iRay][1] << " "
				<< rayPosition[iRay][2] << "\n";

        closestCells = mesh.findHostCellID(rayPosition[iRay], -1);
        std::cout  << "Closest cells = ";

    	for(int i = 0; i < closestCells.size(); i++)
    		std::cout << closestCells[i] << ", ";

    	std::cout  << "\n";
    	std::cout  << "Distance from ray to exit cell = " << distanceToCell << "\n";

    	if(exitCell != -1){
    		std::vector<float> cellPos = mesh.cellCoordinates[exitCell];
    		std::cout << "Next cell (from neighbour) = " << exitCell << "\n"
    				<< "Position of cell (from neighbour search) = "
					<< cellPos[0] << ", " << cellPos[1] << ", " << cellPos[2] << "\n";
    	}
    }


    /*
    if(rayPosition[iRay][0] > mesh.boxSize || rayPosition[iRay][0] < 0 || rayPosition[iRay][1] > mesh.boxSize || rayPosition[iRay][1] < 0 || rayPosition[iRay][2] > mesh.boxSize || rayPosition[iRay][2] < 0){
    	insideDomain[iRay] = false;
    	return -1;
    }
    */


    if(exitCell == -1){
    	warningIssued = true;
    	flagRay[iRay] = true;
    }

    return exitCell;
 }

double Rays::getOvershootDistance(int exitCell, int iRay, double distanceToExit, bool verbose){

	if(insideDomain[iRay] == false)
		return 0.0;

	double distanceRayToExitCellCentre = mesh.getDistanceToCell(rayPosition[iRay], exitCell);
	std::vector<double> positionTmp (3, 0.0);

	double overshoot = distanceRayToExitCellCentre / 10;

	int nIter = 0;
	do {
		for (int i = 0; i < 3; i++)
			positionTmp[i] = rayPosition[iRay][i] + rayDirection[iRay][i] * (distanceToExit + overshoot);

		nIter +=1;
		if(nIter == maxnIter)
			flagRay[iRay] = true;

		overshoot /= 2;

	} while ((mesh.findHostCellID(positionTmp, exitCell)[0] != exitCell) && (nIter < maxnIter));  // todo maybe replace by a check of the entire list, but a priori there should only be one cell as we no longer are on an edge

    if(verbose)
    	std::cout << "Ray position (before update) = " << rayPosition[iRay][0] << ", "  << rayPosition[iRay][1] << ", " << rayPosition[iRay][2] << "\n";

    return overshoot;
}


bool Rays::updateRayAndIsMaxReached(int iCell, int iRay, double& distanceToExit){

	int filter = 1;
	if(flowFilter != 0)
		filter = getFilterForVelocity(flowFilter, iCell, iRay);

	double newColumnDensity  = columnDensity[iRay]  + distanceToExit * mesh.cellDensity[iCell] * filter;
	double newColumnVelocity = columnVelocity[iRay] + distanceToExit * mesh.cellDensity[iCell] * filter * (rayDirection[iRay][0] * mesh.cellVelocities[iCell][0] + rayDirection[iRay][1] * mesh.cellVelocities[iCell][1] + rayDirection[iRay][2] * mesh.cellVelocities[iCell][2]);
	double newDistanceTravelled = distanceTravelled[iRay] + distanceToExit;

	if((newColumnDensity >= maxColumn) && (distanceToExit > 0.0)){

		double fractionalColumn = (maxColumn - columnDensity[iRay])/(newColumnDensity - columnDensity[iRay]);
		distanceToExit *= fractionalColumn;

		columnDensity[iRay]     += distanceToExit * mesh.cellDensity[iCell] * filter;
		columnVelocity[iRay]    += distanceToExit * mesh.cellDensity[iCell] * filter * (rayDirection[iRay][0] * mesh.cellVelocities[iCell][0] + rayDirection[iRay][1] * mesh.cellVelocities[iCell][1] + rayDirection[iRay][2] * mesh.cellVelocities[iCell][2]);
		distanceTravelled[iRay] += distanceToExit;

		return true;
	}
	else if((newDistanceTravelled >= maxRadius) && (distanceToExit > 0.0)){

		double fractionalDistance = (maxRadius - distanceTravelled[iRay])/(newDistanceTravelled - distanceTravelled[iRay]);
		distanceToExit *= fractionalDistance;

		columnDensity[iRay]     += distanceToExit * mesh.cellDensity[iCell] * filter;
		columnVelocity[iRay]    += distanceToExit * mesh.cellDensity[iCell] * filter * (rayDirection[iRay][0] * mesh.cellVelocities[iCell][0] + rayDirection[iRay][1] * mesh.cellVelocities[iCell][1] + rayDirection[iRay][2] * mesh.cellVelocities[iCell][2]);
		distanceTravelled[iRay] += distanceToExit;
		return true;

	}

	columnDensity[iRay]     = newColumnDensity;
	columnVelocity[iRay]    = newColumnVelocity;
	distanceTravelled[iRay] = newDistanceTravelled;

	return false;
}


bool Rays::shouldRayBeTerminated(int iRay, double distanceToExit){

	 if(distanceToExit > mesh.boxSize){
		std::cout << "Distance to exit larger than box size; ray stopped!" << "\n";
		return true;
	 }

	 if (distanceToExit < 1e-10)
	 	std::cout << "distanceToExit = " << distanceToExit << "You seem to have ended up on an edge; how did you do that?!" << "\n";

	 return false;
}

int Rays::modifyExitCellIfOnInterface(int iCell, int iRay, int exitCell, double distanceToExit, bool verbose){

	double distanceToExitTmp;
	std::vector<float> cellPos = mesh.cellCoordinates[iCell];
	std::vector<double> normalVector (3, 0.0);
	std::vector<double> pointOnInterface (3, 0.0);
	std::vector<double> positionTmp (3, 0.0);

	for (int i = 0; i < 3; i++)
		 positionTmp[i] = rayPosition[iRay][i] + rayDirection[iRay][i] * distanceToExit;


	 int targetCell = mesh.findHostCellID(positionTmp, exitCell)[0];

	 if(targetCell != exitCell && targetCell != iCell){

		 flagRay[iRay] = true;
		 if(verbose)
			 std::cout << "Here there's an issue. iCell = " << iCell << ", exitCell = " << exitCell << ", wants to be in cell = " <<  targetCell << std::endl;

		 for (int neighbour : mesh.neighbourList[iCell]){
			 std::vector<float> neighbourPos = mesh.cellCoordinates[neighbour];

			 for (int i = 0; i < 3; i++){
				 normalVector[i] = neighbourPos[i] - cellPos[i];
				 pointOnInterface[i] = 0.5 * (cellPos[i] + neighbourPos[i]);
			 }

			 double denominator = normalVector[0] * rayDirection[iRay][0] + normalVector[1] * rayDirection[iRay][1] + normalVector[2] * rayDirection[iRay][2];

			 if(denominator <= 0.)
				 continue;

			 distanceToExitTmp = normalVector[0] * (pointOnInterface[0] - rayPosition[iRay][0])
					 + normalVector[1] * (pointOnInterface[1] - rayPosition[iRay][1])
					 + normalVector[2] * (pointOnInterface[2] - rayPosition[iRay][2]);

			 distanceToExitTmp = distanceToExitTmp/denominator;

			 if(verbose)
				 std::cout << "cell index = " << neighbour << ", distance = " << distanceToExitTmp << std::endl;

			 if((distanceToExitTmp < 0) && (distanceToExitTmp > -minTolerance))
				 exitCell = neighbour;
			 else if ((neighbour == targetCell) && (fabs(distanceToExitTmp - distanceToExit) < minTolerance))
				 exitCell = neighbour;

			 	 }
	 }

	 return exitCell;
}


    mesh.cellFlux[exitCell] = exp(-columnDensity[iRay] * 1000);

void Rays::findExitCellAndSetDistance(int iCell, int iRay, int& exitCell, double& distanceToExit, bool verbose){
	 double distanceToExitTmp;
	 std::vector<float> cellPos = mesh.cellCoordinates[iCell];
	 std::vector<double> normalVector (3, 0.0);
	 std::vector<double> pointOnInterface (3, 0.0);

	 for (int neighbour : mesh.neighbourList[iCell]){

		 std::vector<float> neighbourPos = mesh.cellCoordinates[neighbour];

		 if(verbose)
			 std::cout << "Neighbour Index = " << neighbour << ", Neighbour Position = " << neighbourPos[0] << ", " << neighbourPos[1] << ", " << neighbourPos[2] << "\n";

		 for (int i = 0; i < 3; i++){
			 normalVector[i] = neighbourPos[i] - cellPos[i];
			 pointOnInterface[i] = 0.5 * (cellPos[i] + neighbourPos[i]);
		 }

		 double denominator = normalVector[0] * rayDirection[iRay][0] + normalVector[1] * rayDirection[iRay][1] + normalVector[2] * rayDirection[iRay][2];
		 if(denominator <= 0.)
			 continue;

		 distanceToExitTmp = normalVector[0] * (pointOnInterface[0] - rayPosition[iRay][0])
				 + normalVector[1] * (pointOnInterface[1] - rayPosition[iRay][1])
				 + normalVector[2] * (pointOnInterface[2] - rayPosition[iRay][2]);
		 distanceToExitTmp = distanceToExitTmp/denominator;


		 if(verbose)
			 std::cout << std::scientific << std::setprecision(12) << "Distance to Neighbour = " << distanceToExitTmp << " denominator = " << denominator << "\n";


		 if((distanceToExitTmp <= distanceToExit) && (distanceToExitTmp > 0)){
			 distanceToExit = distanceToExitTmp;
			 exitCell = neighbour;

			 if(verbose)
				 std::cout << "New Neighbour candidate = " << neighbour << "\n";


		 }
	 }
}


void Rays::outputResults(std::string& ofileName) {
    std::ofstream outputFile(ofileName);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening file for output!" << std::endl;
        return;
    }

    // Print column headers with fixed widths
    outputFile << std::left
               << std::setw(6)  << "Ray"
               << std::setw(10) << "Theta"
               << std::setw(10) << "Phi"
               << std::setw(15) << "Column"
               << std::setw(15) << "Distance"
               << std::setw(6)  << "Flag"
               << std::setw(10) << "NumCells"
               << std::setw(10) << "LastVisit"
               << std::endl;

    for (int i = 0; i < nRays; ++i) {
        outputFile << std::left
                   << std::setw(6)  << i
                   << std::setw(10) << std::fixed << std::setprecision(3) << theta[i]
                   << std::setw(10) << std::fixed << std::setprecision(3) << phi[i]
                   << std::setw(15) << std::fixed << std::setprecision(6) << columnDensity[i]
                   << std::setw(15) << std::fixed << std::setprecision(6) << distanceTravelled[i]
                   << std::setw(6)  << flagRay[i]
                   << std::setw(10) << numTraversedCells[i]
				   << std::setw(10) << lastVisitedCell[i]
                   << std::endl;
    }

    outputFile.close();
    std::cout << "Results have been written to '" << ofileName << "'" << std::endl;
}



void Rays::doRayTracing(){

	std::vector<double> oldRayPosition (3, 0.0);
	int iCellOld = -1;
	bool debug = false;

	for(int iRay = 0; iRay < nRays; iRay++){
		std::cout << "ray " << iRay << std::endl;
		int iCell = startCell;

		while(insideDomain[iRay]){

			iCellOld = iCell;
			for(int i = 0; i < 3; i++)
				oldRayPosition[i] = rayPosition[iRay][i];

			warningIssued = false;

			visitedCells[iRay].push_back(iCell);
			numTraversedCells[iRay] += 1;

			iCell = travelToNextCell(iCellOld, iRay, false);

			if(debug){
				if(warningIssued){

					for(int i = 0; i < 3; i++)
						rayPosition[iRay][i] = oldRayPosition[i];

					iCell = travelToNextCell(iCellOld, iRay, true);
					insideDomain[iRay] = false;
				}
			}
		}
	}
}


Rays::~Rays() {
	// TODO Auto-generated destructor stub
}
