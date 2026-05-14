/*
 * Rays.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Rays.h"
#include "Mesh.h"

Rays::Rays(double ionisationCrossSectionHI, double ionisationCrossSectionHeI, double ionisationCrossSectionHeII,
		double maxRadius, int64_t Nside, Mesh& mesh) :
		ionisationCrossSectionHI(ionisationCrossSectionHI), ionisationCrossSectionHeI(ionisationCrossSectionHeI), ionisationCrossSectionHeII(ionisationCrossSectionHeII), maxRadius(maxRadius), Nside(Nside), nRays(12 * Nside * Nside), mesh(mesh) {

	ionisationCrossSectionHI_inInternalUnits   = ionisationCrossSectionHI   / mesh.protonMass * mesh.unitMass / mesh.unitLength / mesh.unitLength;
	ionisationCrossSectionHeI_inInternalUnits  = ionisationCrossSectionHeI  / 4.0 / mesh.protonMass * mesh.unitMass / mesh.unitLength / mesh.unitLength;
	ionisationCrossSectionHeII_inInternalUnits = ionisationCrossSectionHeII / 4.0 / mesh.protonMass * mesh.unitMass / mesh.unitLength / mesh.unitLength;
	dustAbsorptionOpacity_inInternalUnits      = dustAbsorptionOpacity  * mesh.unitMass / mesh.unitLength / mesh.unitLength;

	rayDirection = std::vector<std::array<double,3>>(nRays);
	rayPosition  = std::vector<std::array<double,3>>(nRays);

	theta     = std::vector<double>(nRays, 0.0);
	phi       = std::vector<double>(nRays, 0.0);
	rayWeight = std::vector<double>(nRays, 1.0 / nRays);

	columnHI          = std::vector<double>(nRays, 0.0);
	columnDust        = std::vector<double>(nRays, 0.0);
	distanceTravelled = std::vector<double>(nRays, 0.0);
	finalLuminosity   = std::vector<double>(nRays, 0.0);
	insideDomain      = std::vector<bool>(nRays, true);
	flagRay           = std::vector<bool>(nRays, false);

	visitedCellColumn   = std::vector<std::vector<double>>(nRays);
	visitedCellDistance = std::vector<std::vector<double>>(nRays);
	visitedCells        = std::vector<std::vector<int>>(nRays);

	initializeDirections();
}

void Rays::initializeDirections() {

    rayDirection.resize(nRays);
    theta.resize(nRays);
    phi.resize(nRays);
    rayWeight.resize(nRays, 1.0 / nRays);

    for(int iRay = 0; iRay < nRays; ++iRay){

        double thetaTmp, phiTmp;

        pix2ang_ring(Nside, iRay, &thetaTmp, &phiTmp);

        theta[iRay] = thetaTmp;
        phi[iRay]   = phiTmp;

        rayDirection[iRay][0] = std::sin(thetaTmp) * std::cos(phiTmp);
        rayDirection[iRay][1] = std::sin(thetaTmp) * std::sin(phiTmp);
        rayDirection[iRay][2] = std::cos(thetaTmp);
    }
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

	 if(exitCell == -1)
		 insideDomain[iRay] = false;

	 if(insideDomain[iRay]){

		 if(updateRayAndIsMaxReached(iCell, iRay, distanceToExit))
			 insideDomain[iRay] = false;

		 visitedCellColumn[iRay].push_back(distanceToExit * mesh.getDensity(iCell));
		 visitedCells[iRay].push_back(iCell);

		 mesh.cellVisitsByRay[iCell] += 1;

		 overshoot = getOvershootDistance(exitCell, iRay, distanceToExit, verbose);

		 if(updateRayAndIsMaxReached(exitCell, iRay, overshoot))
			 insideDomain[iRay] = false;

		 if (!visitedCellColumn[iRay].empty()){
			 visitedCellColumn[iRay].back() += overshoot * mesh.getDensity(exitCell);
		 }
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
	std::array<double,3> positionTmp{0.0, 0.0, 0.0};

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
		visitedCellDistance[iRay].back() += distanceToExit;

		return true;
	}

	distanceTravelled[iRay] = newDistanceTravelled;
	visitedCellDistance[iRay].push_back(newDistanceTravelled);

	return false;
}


int Rays::modifyExitCellIfOnInterface(int iCell, int iRay, int exitCell, double distanceToExit, bool verbose){

	double distanceToExitTmp;
	std::vector<float> cellPos = mesh.cellCoordinates[iCell];
	std::vector<double> normalVector (3, 0.0);
	std::vector<double> pointOnInterface (3, 0.0);
	std::array<double,3> positionTmp{};

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
    std::cout << "Ray properties have been written to '" << ofileName << "'" << std::endl;
}

double Rays::distanceSquared(std::vector<float>& a, std::vector<float>& b){
	return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
}

void Rays::updateColumnAndFlux(int iRay){

	double NdotFinal  = rayWeight[iRay];

	for (int i = 0; i < visitedCells[iRay].size(); i++){
		int iCell        = visitedCells[iRay][i];

		double xHII = mesh.xH_old[iCell];

		xHII = std::max(0.0, std::min(1.0, xHII));
		double neutral = 1.0 - xHII;

		double yHe       = mesh.getHeIIFraction(iCell);
		double zHe       = mesh.getHeIIIFraction(iCell);
		double neutralHe = 1.0 - yHe - zHe;

		double dColumnHI     = visitedCellColumn[iRay][i] * mesh.xHydrogen * neutral;
		double dColumnHeI    = visitedCellColumn[iRay][i] * mesh.yHelium   * neutralHe;
		double dColumnHeII   = visitedCellColumn[iRay][i] * mesh.yHelium   * yHe;
		double dColumnDust   = visitedCellColumn[iRay][i]; //todo

		const double tauHI   = ionisationCrossSectionHI_inInternalUnits    * dColumnHI;
		const double tauHeI  = ionisationCrossSectionHeI_inInternalUnits   * dColumnHeI;
		const double tauHeII = ionisationCrossSectionHeII_inInternalUnits  * dColumnHeII;
		const double tauDust = dustAbsorptionOpacity_inInternalUnits       * dColumnDust;

		const double dtau = 0.;//tauHI + tauHeI + tauHeII + tauDust;
		const double fabs = 1.0 - std::exp(-dtau);

		mesh.cellIncomingPhotonRate[iCell] += NdotFinal;

		columnHI[iRay]   += dColumnHI;
		columnDust[iRay] += dColumnDust;

		NdotFinal *= 1.0 - fabs;
	}

	finalLuminosity[iRay] = NdotFinal;

}

void Rays::calculateRays(){

    for(int iCell = 0; iCell < mesh.numCells; ++iCell){

        for(int iRay = 0; iRay < nRays; ++iRay){

            resetRay(iRay, iCell);

            int iCellCurrent = iCell;

            while(insideDomain[iRay])
                iCellCurrent = travelToNextCell(iCellCurrent, iRay, false);

            updateColumnAndFlux(iRay);
        }
    }
}

void Rays::resetRay(int iRay, int startCell){

    insideDomain[iRay]      = true;
    flagRay[iRay]           = false;

    columnHI[iRay]          = 0.0;
    columnDust[iRay]        = 0.0;
    distanceTravelled[iRay] = 0.0;
    finalLuminosity[iRay]   = 0.0;

    visitedCellColumn[iRay].clear();
    visitedCellDistance[iRay].clear();
    visitedCells[iRay].clear();

    for(int i = 0; i < 3; ++i)
        rayPosition[iRay][i] = mesh.cellCoordinates[startCell][i];
}

Rays::~Rays() {
}
