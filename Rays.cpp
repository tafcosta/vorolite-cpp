/*
 * Rays.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Rays.h"
#include "Mesh.h"

Rays::Rays(double ionisationCrossSection, double maxRadius, std::vector<double> sourcePosition, double lumTotal, Mesh& mesh) : ionisationCrossSection(ionisationCrossSection), maxRadius(maxRadius), sourcePosition(sourcePosition), lumTotal(lumTotal), mesh(mesh) {

	ionisationCrossSection_inInternalUnits = ionisationCrossSection / mesh.protonMass * mesh.unitMass / mesh.unitLength / mesh.unitLength;

	startCell = mesh.findHostCellID(sourcePosition, -1)[0];

	setNumRays();

	rayTargetCell = std::vector<int> (nRays);
    for (int iRay = 0; iRay < nRays; ++iRay)
    	rayTargetCell[iRay] = iRay;

	rayDirection = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	theta        = std::vector<double>(nRays, 0.0);
	phi          = std::vector<double>(nRays, 0.0);
	rayWeight    = std::vector<double>(nRays, 0.0);
	initializeDirections();
	assignToHealpix();

	rayPosition = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	initializePositions();

	columnHI          = std::vector<double>(nRays, 0.0);
	distanceTravelled = std::vector<double>(nRays, 0.0);
	insideDomain      = std::vector<bool>(nRays, true);
	flagRay           = std::vector<bool>(nRays, false);

	visitedCellColumn = std::vector<std::vector<double>>(nRays);
	visitedCells      = std::vector<std::vector<int>>(nRays);

}

double Rays::getLuminosity(double time){
	//return (time < 0.01? 0 : (time > 0.02 ? 0. : lumTotal));]
	return lumTotal;
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

    	if(rayTargetCell[iRay] != startCell){
    		xDistance = cellPos[0] - mesh.cellCoordinates[startCell][0];
    		yDistance = cellPos[1] - mesh.cellCoordinates[startCell][1];
    		zDistance = cellPos[2] - mesh.cellCoordinates[startCell][2];
    	}
    	else
    	{
    		xDistance = cellPos[0] - sourcePosition[0];
    		yDistance = cellPos[1] - sourcePosition[1];
    		zDistance = cellPos[2] - sourcePosition[2];
    	}

    	rDistance = std::sqrt(xDistance*xDistance + yDistance*yDistance + zDistance*zDistance);

    	if(rDistance > 0.){
    		rayDirection[iRay][0] = xDistance / rDistance;
    		rayDirection[iRay][1] = yDistance / rDistance;
    		rayDirection[iRay][2] = zDistance / rDistance;

    		phi[iRay]   = std::atan2(yDistance, xDistance);
    		theta[iRay] = std::acos(zDistance / rDistance);
    	}
    }
}

void Rays::assignToHealpix() {
	int64_t healpixNside = 32;
	int64_t nPix = nside2npix(healpixNside);
	std::vector<int> raysPerPixel(nPix, 0);

	std::vector<int64_t> rayToPixel(nRays);
	for (int i = 0; i < nRays; ++i) {
		long iPix;
		ang2pix_ring(healpixNside, theta[i], phi[i], &iPix);

		rayToPixel[i] = static_cast<int64_t>(iPix);
		raysPerPixel[iPix]++;
	}

	const double omegaPix = 4.0 * M_PI / static_cast<double>(nPix);

	for (int i = 0; i < nRays; ++i) {
		int64_t iPix = rayToPixel[i];
		rayWeight[i] = 1.0 / raysPerPixel[iPix] * omegaPix / (4.0 * M_PI) ;
	}
}


void Rays::initializePositions() {
    for (int iRay = 0; iRay < nRays; ++iRay)
        for (int i = 0; i < 3; ++i)
        	rayPosition[iRay][i] = 	mesh.cellCoordinates[startCell][i];

    std::cout << "Source Position changed to = " << mesh.cellCoordinates[startCell][0] << " " << mesh.cellCoordinates[startCell][1] << " " << mesh.cellCoordinates[startCell][2] << std::endl;
}


 int Rays::travelToNextCell(int iCell, int iRay, bool verbose){
	 double distanceToExit = std::numeric_limits<double>::max();
	 double overshoot = 0.;
	 int exitCell     = -1;

	 if(verbose){
		 std::vector<float> cellPos = mesh.cellCoordinates[iCell];
		 double distanceBetRayAndCell = sqrt((cellPos[0] - rayPosition[iRay][0]) * (cellPos[0] - rayPosition[iRay][0]) + (cellPos[1] - rayPosition[iRay][1]) * (cellPos[1] - rayPosition[iRay][1]) + (cellPos[2] - rayPosition[iRay][2]) * (cellPos[2] - rayPosition[iRay][2]));
	 
		 std::cout << "Host Index = " << iCell << " Host Cell Position = " << cellPos[0] << ", " << cellPos[1] << ", " << cellPos[2] << " Distance from Ray to Cell (before update) = " << distanceBetRayAndCell << "\n";
	 }

	 exitCell =  findExitCellAndSetDistance(iCell, iRay, exitCell, distanceToExit, verbose);
	 exitCell = modifyExitCellIfOnInterface(iCell, iRay, exitCell, distanceToExit, verbose);

	 if (distanceToExit < 1e-10)
	     std::cerr << "Warning: distanceToExit = " << distanceToExit
	               << ". You seem to have ended up on an edge; how did you do that?!" << std::endl;

	 if(exitCell == -1){
	     std::cerr << "No exit cell found. You need a larger domain buffer size." << std::endl;
		 insideDomain[iRay] = false;
	 }

	 if(insideDomain[iRay]){

		 if(updateRayAndIsMaxReached(iCell, iRay, distanceToExit))
			 insideDomain[iRay] = false;

		 visitedCellColumn[iRay].push_back(distanceToExit * mesh.getDensity(iCell));
		 visitedCells[iRay].push_back(iCell);

		 overshoot = getOvershootDistance(exitCell, iRay, distanceToExit, verbose);

		 if(updateRayAndIsMaxReached(exitCell, iRay, overshoot))
			 insideDomain[iRay] = false;

		 if (!visitedCellColumn[iRay].empty())
			 visitedCellColumn[iRay].back() += overshoot * mesh.getDensity(exitCell);
	 }


	 updateRayPosition(iRay, distanceToExit + overshoot);


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


    if(exitCell == -1){
    	warningIssued = true;
    	flagRay[iRay] = true;
    }

	 if(insideDomain[iRay] == false)
		 return -1;


    return exitCell;
 }

 void Rays::updateRayPosition(int iRay, double distance){
	 for (int i = 0; i < 3; i++)
		 rayPosition[iRay][i] += rayDirection[iRay][i] * distance;
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

	double newDistanceTravelled = distanceTravelled[iRay] + distanceToExit;

	if(distanceTravelled[iRay] + distanceToExit > maxRadius){

		double fractionalDistance = (newDistanceTravelled - distanceTravelled[iRay])/(newDistanceTravelled - maxRadius);
		distanceToExit *= fractionalDistance;

		distanceTravelled[iRay] += distanceToExit;
		return true;
	}

	distanceTravelled[iRay] = newDistanceTravelled;
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


int Rays::findExitCellAndSetDistance(int iCell, int iRay, int& exitCell, double& distanceToExit, bool verbose){
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

	 return exitCell;
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
               << std::setw(10) << "LastVisit"
               << std::endl;

    for (int i = 0; i < nRays; ++i) {
        outputFile << std::left
                   << std::setw(6)  << i
                   << std::setw(10) << std::fixed << std::setprecision(3) << theta[i]
                   << std::setw(10) << std::fixed << std::setprecision(3) << phi[i]
                   << std::setw(15) << std::fixed << std::setprecision(6) << columnHI[i]
                   << std::setw(15) << std::fixed << std::setprecision(6) << distanceTravelled[i]
                   << std::setw(6)  << flagRay[i]
				   << std::setw(10) << rayTargetCell[i]
                   << std::endl;
    }

    outputFile.close();
    std::cout << "Results have been written to '" << ofileName << "'" << std::endl;
}

double Rays::distanceSquared(std::vector<float>& a, std::vector<float>& b){
	return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
}

void Rays::updateColumnAndFlux(int iRay, double time, double dtime){
	columnHI[iRay] = 0.;

	if(timeDependent){

		int j = 0;
		double columnHIIindTime = 0.;
		std::vector<double> fluxOfRay (visitedCells[iRay].size(), 0.0);

		for (int i = 0; i < visitedCells[iRay].size(); i++){

			if(i > 0)
				columnHIIindTime += visitedCellColumn[iRay][i] * (1 - mesh.getHIIFraction(visitedCells[iRay][i])) / 2.0
					+ visitedCellColumn[iRay][i-1] * (1 - mesh.getHIIFraction(visitedCells[iRay][i-1])) / 2.0;

			//todo: should use linear interpolation instead
			while((distanceSquared(mesh.cellCoordinates[visitedCells[iRay][j]], mesh.cellCoordinates[visitedCells[iRay][i]]) > speedOfLightInternal * speedOfLightInternal * dtime * dtime)
					&& (j < visitedCells[iRay].size() - 2) && (j < i)){

				columnHIIindTime -= visitedCellColumn[iRay][j] * (1 - mesh.getHIIFraction(visitedCells[iRay][j])) / 2.0;
						if(j < visitedCells[iRay].size() - 1)
							columnHIIindTime -= visitedCellColumn[iRay][j+1] * (1 - mesh.getHIIFraction(visitedCells[iRay][j+1])) / 2.0;

				j++;
			}

			fluxOfRay[i] = (j == 0? getLuminosity(time - dtime/2.0) * rayWeight[iRay]: mesh.getFluxOfRayInCell(iRay, j)) * std::exp(-ionisationCrossSection_inInternalUnits * columnHIIindTime);
			mesh.cellFlux[visitedCells[iRay][i]] += fluxOfRay[i];

			if(visitedCells[iRay][i] == rayTargetCell[iRay])
				mesh.cellLocalColumn[rayTargetCell[iRay]] = visitedCellColumn[iRay][i];
		}

		for (int i = 0; i < visitedCells[iRay].size(); i++)
			mesh.setFluxOfRayInCell(iRay, i, fluxOfRay[i]);

	}
	else
	{

		for (int i = 0; i < visitedCells[iRay].size(); i++){
			mesh.cellIncomingFlux[visitedCells[iRay][i]] += getLuminosity(0.0) * rayWeight[iRay] * std::exp(-ionisationCrossSection_inInternalUnits * columnHI[iRay]) ;

			columnHI[iRay] += visitedCellColumn[iRay][i] * (1 - mesh.getHIIFraction(visitedCells[iRay][i]));
		    mesh.cellFlux[visitedCells[iRay][i]] += getLuminosity(0.0) * rayWeight[iRay] * std::exp(-ionisationCrossSection_inInternalUnits * columnHI[iRay]);

		    if(visitedCells[iRay][i] == rayTargetCell[iRay])
				mesh.cellLocalColumn[rayTargetCell[iRay]] = visitedCellColumn[iRay][i];
		}
	}

}

void Rays::printVisitedCellColumn(){

	for(int iRay = 0; iRay < nRays; iRay++)
		for (int i = 0; i < visitedCells[iRay].size(); i++)
			if(iRay == 1106)
				std::cout << visitedCells[iRay][i] << " " << visitedCellColumn[iRay][i]  << std::endl;
}

void Rays::calculateRays(){

	for(int iRay = 0; iRay < nRays; iRay++){
		int iCell = startCell;
		while(insideDomain[iRay])
			iCell = travelToNextCell(iCell, iRay, false);

		mesh.resizeFluxOfRayInCell(iRay, visitedCells[iRay].size());
		}
}


void Rays::doRadiativeTransfer(double time, double dtime){
	for(int iRay = 0; iRay < nRays; iRay++)
		updateColumnAndFlux(iRay, time, dtime);
}


Rays::~Rays() {
	// TODO Auto-generated destructor stub
}
