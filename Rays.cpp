/*
 * Rays.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Rays.h"
#include "Mesh.h"

Rays::Rays(int nRays, double maxRadius, std::vector<double> sourcePosition, int flowFilter, double maxColumn, Mesh& mesh) : numRays(nRays), maxRadius(maxRadius), sourcePosition(sourcePosition), flowFilter(flowFilter), maxColumn(maxColumn), mesh(mesh) {
	startCell = mesh.findHostCellID(sourcePosition, -1)[0];

	rayDirection = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	theta = std::vector<double>(nRays, 0.0);
	phi = std::vector<double>(nRays, 0.0);
	initializeDirections();

	rayPosition = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	initializePositions();

	visitedCells      = std::vector<std::vector<int>>(nRays);
	columnDensity     = std::vector<double>(nRays, 0.0);
	columnVelocity    = std::vector<double>(nRays, 0.0);

	distanceTravelled = std::vector<double>(nRays, 0.0);
	insideDomain      = std::vector<bool>(nRays, true);
	flagRay           = std::vector<bool>(nRays, false);
	numTraversedCells = std::vector<int>(nRays, 0.0);

	warningIssued = false;
}

void Rays::initializeDirections() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int iRay = 0; iRay < numRays; ++iRay) {
        phi[iRay] = dis(gen) * 2.0 * M_PI;
        theta[iRay] = std::acos(1.0 - 2.0 * dis(gen));

        rayDirection[iRay][0] = std::cos(phi[iRay]) * std::sin(theta[iRay]);
        rayDirection[iRay][1] = std::sin(phi[iRay]) * std::sin(theta[iRay]);
        rayDirection[iRay][2] = std::cos(theta[iRay]);
    }
}

void Rays::initializePositions() {
    for (int iRay = 0; iRay < numRays; ++iRay)
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


	 if(distanceToExit > mesh.boxSize){
		std::cout << "BOO!" << " " << exitCell << " " << distanceToExit << std::endl;
	 }


	 exitCell = modifyExitCellIfOnInterface(iCell, iRay, exitCell, distanceToExit, verbose);

	 if(shouldRayBeTerminated(iRay, distanceToExit))
		 insideDomain[iRay] = false;

	 if(insideDomain[iRay]){

		 if(updateColumnAndIsMaxReached(iCell, iRay, distanceToExit))
			 insideDomain[iRay] = false;

		 overshoot = getOvershootDistance(exitCell, iRay, distanceToExit, verbose);

		 if(updateColumnAndIsMaxReached(exitCell, iRay, overshoot))
			 insideDomain[iRay] = false;
	 }


	 distanceTravelled[iRay] += distanceToExit + overshoot;

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







    if(rayPosition[iRay][0] > mesh.boxSize || rayPosition[iRay][0] < 0 || rayPosition[iRay][1] > mesh.boxSize || rayPosition[iRay][1] < 0 || rayPosition[iRay][2] > mesh.boxSize || rayPosition[iRay][2] < 0 || distanceTravelled[iRay] >= maxRadius || columnDensity[iRay] >= maxColumn){
    	insideDomain[iRay] = false;
    	return -1;
    }

    if(exitCell == -1){
    	warningIssued = true;
    	flagRay[iRay] = true;
    }

    return exitCell;
 }

double Rays::getOvershootDistance(int exitCell, int iRay, double distanceToExit, bool verbose){

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


bool Rays::updateColumnAndIsMaxReached(int iCell, int iRay, double& distanceToExit){

	int filter = 1;
	if(flowFilter != 0)
		filter = getFilterForVelocity(flowFilter, iCell, iRay);

	double newColumnDensity  = columnDensity[iRay]  + distanceToExit * mesh.cellDensity[iCell] * filter;
	double newColumnVelocity = columnVelocity[iRay] + distanceToExit * mesh.cellDensity[iCell] * filter * (rayDirection[iRay][0] * mesh.cellVelocities[iCell][0] + rayDirection[iRay][1] * mesh.cellVelocities[iCell][1] + rayDirection[iRay][2] * mesh.cellVelocities[iCell][2]);

	if((newColumnDensity >= maxColumn) && (maxColumn > 0.)){

		double fractionalColumn = (maxColumn - columnDensity[iRay])/(newColumnDensity - columnDensity[iRay]);

		distanceToExit *= fractionalColumn;

		columnDensity[iRay]     += distanceToExit * mesh.cellDensity[iCell] * filter;
		columnVelocity[iRay]    += distanceToExit * mesh.cellDensity[iCell] * filter * (rayDirection[iRay][0] * mesh.cellVelocities[iCell][0] + rayDirection[iRay][1] * mesh.cellVelocities[iCell][1] + rayDirection[iRay][2] * mesh.cellVelocities[iCell][2]);
		return true;

	} else {
		columnDensity[iRay] = newColumnDensity;
		columnVelocity[iRay] = newColumnVelocity;
	}

	return false;

}


bool Rays::shouldRayBeTerminated(int iRay, double distanceToExit){

	 if(distanceToExit > mesh.boxSize){
		std::cout << "Distance to exit larger than box size; stopped!" << "\n";
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
               << std::setw(15) << "Velocity"
               << std::setw(15) << "Distance"
               << std::setw(6)  << "Flag"
               << std::setw(10) << "NumCells"
               << std::endl;

    for (int i = 0; i < numRays; ++i) {
        outputFile << std::left
                   << std::setw(6)  << i
                   << std::setw(10) << std::fixed << std::setprecision(3) << theta[i]
                   << std::setw(10) << std::fixed << std::setprecision(3) << phi[i]
                   << std::setw(15) << std::fixed << std::setprecision(6) << columnDensity[i]
                   << std::setw(15) << std::fixed << std::setprecision(6) << columnVelocity[i]
                   << std::setw(15) << std::fixed << std::setprecision(6) << distanceTravelled[i]
                   << std::setw(6)  << flagRay[i]
                   << std::setw(10) << numTraversedCells[i]
                   << std::endl;
    }

    outputFile.close();
    std::cout << "Results have been written to '" << ofileName << "'" << std::endl;
}


 double Rays::getFilterForVelocity(int flowFilter, int cellIndex, int iRay){

	 const auto& rayDir = rayDirection[iRay];
	 const auto& velocity = mesh.cellVelocities[cellIndex];

	 double dotProduct = rayDir[0] * velocity[0] + rayDir[1] * velocity[1] + rayDir[2] * velocity[2];

	 if(flowFilter != 0)
		 return (flowFilter * dotProduct > 0.0) ? 1.0 : 0.0;

	 return 1.0;
 }


void Rays::doRayTracing(){

	std::vector<double> oldRayPosition (3, 0.0);
	int iCellOld = -1;
	bool debug = false;

	for(int iRay = 0; iRay < numRays; iRay++){
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

		if(columnDensity[iRay] > 0.)
			columnVelocity[iRay] = columnVelocity[iRay] / columnDensity[iRay];
		else
			columnVelocity[iRay] = 0.;
	}
}


Rays::~Rays() {
	// TODO Auto-generated destructor stub
}
