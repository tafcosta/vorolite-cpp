/*
 * Rays.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Rays.h"
#include "Mesh.h"

Rays::Rays(int nRays, std::vector<double> sourcePosition, Mesh& mesh) : numRays(nRays), sourcePosition(sourcePosition), mesh(mesh) {
	rayDirection = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	initializeDirections();
	rayPosition = std::vector<std::vector<double>>(nRays, sourcePosition);

	opticalDepth = std::vector<double>(nRays, 0.0);
	numTraversedCells = std::vector<double>(nRays, 0.0);
	insideDomain = std::vector<bool>(nRays, true);

	startCell = mesh.findHostCellID(sourcePosition);

}

void Rays::initializeDirections() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int iRay = 0; iRay < numRays; ++iRay) {
        double phi = dis(gen) * 2.0 * M_PI;
        double theta = std::acos(1.0 - 2.0 * dis(gen));

        rayDirection[iRay][0] = 0;//std::cos(phi) * std::sin(theta);
        rayDirection[iRay][1] = 0;//std::sin(phi) * std::sin(theta);
        rayDirection[iRay][2] = 1;//std::cos(theta);
    }
}

 int Rays::findNextCell(int iCell, int iRay){

	 double distanceToExit = std::numeric_limits<double>::max();
	 double distanceToExitTmp;
	 int exitCell  = -1;

	 std::vector<float> cellPos = mesh.cellCoordinates[iCell];

	 std::cout << "Ray Pos = " << rayPosition[iRay][0] << " " << rayPosition[iRay][1] << " " << rayPosition[iRay][2] << std::endl;
	 std::cout << "number of neighbours = " << mesh.neighbourList[iCell].size() << std::endl;

	 for (int neighbour : mesh.neighbourList[iCell]){
		 std::vector<float> neighbourPos = mesh.cellCoordinates[neighbour];
		 std::vector<double> normalVector (3, 0.0);
		 std::vector<double> pointOnInterface (3, 0.0);

		 for (int i = 0; i < 3; i++){
		 normalVector[i] = cellPos[i] - neighbourPos[i];
		 pointOnInterface[i]= 0.5 *(cellPos[i] + neighbourPos[i]);
		 }

		 double denominator = normalVector[0] * rayDirection[iRay][0] + normalVector[1] * rayDirection[iRay][1] + normalVector[2] * rayDirection[iRay][2];

		 if(fabs(denominator) > std::numeric_limits<double>::epsilon()){
			 distanceToExitTmp = normalVector[0] * (pointOnInterface[0] - rayPosition[iRay][0])
					 + normalVector[1] * (pointOnInterface[1] - rayPosition[iRay][1])
					 + normalVector[2] * (pointOnInterface[2] - rayPosition[iRay][2]);
	        distanceToExitTmp = distanceToExitTmp/denominator;
		 } else
			 continue;
	        	
		 if((distanceToExitTmp < distanceToExit) & (distanceToExitTmp > std::numeric_limits<double>::epsilon())){
			 distanceToExit = distanceToExitTmp;
			 exitCell = neighbour;
		 }
	 }

  	std::cout <<  distanceToExit << std::endl;

 	if(exitCell == -1){
 		std::cerr << "Error: no exit cell found!" << std::endl;
 		insideDomain[iRay] = false;
 	}

 	for (int i = 0; i < 3; i++)
 		rayPosition[iRay][i] += rayDirection[iRay][i] * distanceToExit * 1.0000001;

 	numTraversedCells[iRay] += 1;

    return exitCell;

 }


void Rays::doRayTracing(){

	for(int iRay = 0; iRay < numRays; iRay++){

		int iCell = startCell;

		while(insideDomain[iRay]){
			iCell = findNextCell(iCell, iRay);
			std::cout << "exit = " << iCell << std::endl;
		}
	}

}


Rays::~Rays() {
	// TODO Auto-generated destructor stub
}

