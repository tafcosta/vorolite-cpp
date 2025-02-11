/*
 * Rays.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Rays.h"
#include "Mesh.h"

Rays::Rays(int nRays, double maxRadius, std::vector<double> sourcePosition, Mesh& mesh) : numRays(nRays), maxRadius(maxRadius), sourcePosition(sourcePosition), mesh(mesh) {
	rayDirection = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	theta = std::vector<double>(nRays, 0.0);
	phi = std::vector<double>(nRays, 0.0);

	initializeDirections();
	rayPosition = std::vector<std::vector<double>>(nRays, sourcePosition);

	columnDensity = std::vector<double>(nRays, 0.0);
	distanceTravelled = std::vector<double>(nRays, 0.0);
	numTraversedCells = std::vector<double>(nRays, 0.0);
	insideDomain = std::vector<bool>(nRays, true);

	startCell = mesh.findHostCellID(sourcePosition, -1);

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

 int Rays::findNextCell(int iCell, int iRay){
	 double distanceToExit = std::numeric_limits<double>::max();
	 double distanceToExitTmp;
	 int exitCell  = -1;
	 bool verbose = false;

	 std::vector<float> cellPos = mesh.cellCoordinates[iCell];

	 if(verbose){
		 std::cout << "Host Index = " << iCell << " Host Cell Position = " << cellPos[0] << ", " << cellPos[1] << ", " << cellPos[2] << std::endl;
		 double distanceBetRayAndCell = sqrt((cellPos[0] - rayPosition[iRay][0]) * (cellPos[0] - rayPosition[iRay][0]) + (cellPos[1] - rayPosition[iRay][1]) * (cellPos[1] - rayPosition[iRay][1]) + (cellPos[2] - rayPosition[iRay][2]) * (cellPos[2] - rayPosition[iRay][2]));
		 std::cout << "Distance from Ray to Cell (before update) = " << distanceBetRayAndCell << std::endl;
	 }

	 for (int neighbour : mesh.neighbourList[iCell]){
		 std::vector<float> neighbourPos = mesh.cellCoordinates[neighbour];

		 if(verbose)
			 std::cout << "Neighbour Index = " << neighbour << ", " << " Neighbour Position = " << neighbourPos[0] << ", " << neighbourPos[1] << ", " << neighbourPos[2] << std::endl;

		 std::vector<double> normalVector (3, 0.0);
		 std::vector<double> pointOnInterface (3, 0.0);

		 for (int i = 0; i < 3; i++){
			 normalVector[i] = neighbourPos[i] - cellPos[i];
			 pointOnInterface[i] = 0.5 * (cellPos[i] + neighbourPos[i]);
		 }


		 double denominator = normalVector[0] * rayDirection[iRay][0] + normalVector[1] * rayDirection[iRay][1] + normalVector[2] * rayDirection[iRay][2];

		 distanceToExitTmp = normalVector[0] * (pointOnInterface[0] - rayPosition[iRay][0])
				 + normalVector[1] * (pointOnInterface[1] - rayPosition[iRay][1])
				 + normalVector[2] * (pointOnInterface[2] - rayPosition[iRay][2]);
		 distanceToExitTmp = distanceToExitTmp/denominator;

		 if(verbose)
			 std::cout << "Distance to Neighbour = " << distanceToExitTmp << std::endl;

		 if(denominator > 0.){

		 } else
			 continue;



		 if((distanceToExitTmp < distanceToExit) && (distanceToExitTmp >= 0.)){


			 if(distanceToExitTmp == 0){

				 if((neighbourPos[0] - rayPosition[iRay][0]) * rayDirection[iRay][0] + (neighbourPos[1] - rayPosition[iRay][1]) * rayDirection[iRay][1] + (neighbourPos[2] - rayPosition[iRay][2]) * rayDirection[iRay][2] < 0)
					 continue;
			 }


			 distanceToExit = distanceToExitTmp;
			 exitCell = neighbour;

			 if(verbose)
				 std::cout << "New Neighbour candidate= " << neighbour << std::endl;

		 }
	 }

	 if(verbose)
		 std::cout << "Exit Cell= " << exitCell << std::endl;

	 if(distanceToExit > mesh.boxSize){
	 	insideDomain[iRay] = false;

	 	distanceToExit = 0.0;
	 	exitCell = -1;
	 }

	 if(verbose){
		 std::cout << "Ray position (before update) = " << rayPosition[iRay][0] << " " << rayPosition[iRay][1] << " " << rayPosition[iRay][2] << std::endl;
		 std::cout << "Cell based on ray position (before update) = " << mesh.findHostCellID(rayPosition[iRay], exitCell) << std::endl;
	 }

	 for (int i = 0; i < 3; i++)
		 rayPosition[iRay][i] += rayDirection[iRay][i] * distanceToExit;

	 columnDensity[iRay] += distanceToExit * mesh.cellDensity[iCell];
	 distanceTravelled[iRay] += distanceToExit;
	 numTraversedCells[iRay] += 1;


	 if(rayPosition[iRay][0] > 1 || rayPosition[iRay][0] < 0 || rayPosition[iRay][1] > 1 || rayPosition[iRay][1] < 0 || rayPosition[iRay][2] > 1 || rayPosition[iRay][2] < 0 || distanceTravelled[iRay] >= maxRadius){
		 insideDomain[iRay] = false;
		 exitCell = -1;
	 } else {
		 if(exitCell == -1)
			 std::cerr << "Warning: inside domain, but no neighbours found. Ray number = " << iRay << std::endl;
	 }


	 if(verbose){
		 if(exitCell != mesh.findHostCellID(rayPosition[iRay], iCell)){
			 std::cout << "Warning: exit cell does not match updated cell hosting ray. Ray position = " << rayPosition[iRay][0] << " " << rayPosition[iRay][1] << " " << rayPosition[iRay][2] << std::endl;
			 std::cout << "Indices = " << exitCell << ", " << mesh.findHostCellID(rayPosition[iRay], iCell) << ", distance = " << distanceToExit << std::endl;
		 }

		 std::cout << "Ray position (after update) = " << rayPosition[iRay][0] << " " << rayPosition[iRay][1] << " " << rayPosition[iRay][2] << std::endl;
		 std::cout << "Column density = " << columnDensity[iRay] << std::endl;

		 std::cout << "Cell based on ray position (after update) = " << mesh.findHostCellID(rayPosition[iRay], iCell) << std::endl;

		 cellPos = mesh.cellCoordinates[exitCell];

		 std::cout << "Next cell from neighbour = " << exitCell << std::endl;
		 std::cout << "Position of cell (from neighbour search) = " << cellPos[0] << ", " << cellPos[1] << ", " << cellPos[2] << std::endl;
		 std::cout << "-----" << std::endl;
	 }

	 return exitCell;

 }

 void Rays::outputResults() {
     std::ofstream outputFile("ray_output.txt");

     if (!outputFile.is_open()) {
         std::cerr << "Error opening file for output!" << std::endl;
         return;
     }

     outputFile << "Ray Number\tTheta (radians)\tPhi (radians)\tColumn Density" << std::endl;

     for (int i = 0; i < numRays; ++i) {
         outputFile << i + 1 << "\t"  // Ray number
                    << std::fixed << std::setprecision(3) << theta[i] << "\t"
                    << std::fixed << std::setprecision(3) << phi[i] << "\t"
                    << std::fixed << std::setprecision(6) << columnDensity[i] << std::endl;
     }

     outputFile.close();

     std::cout << "Results have been written to 'ray_output.txt'" << std::endl;
 }


void Rays::doRayTracing(){

	for(int iRay = 0; iRay < numRays; iRay++){

		int iCell = startCell;

		while(insideDomain[iRay]){
			iCell = findNextCell(iCell, iRay);
		}
	}

}


Rays::~Rays() {
	// TODO Auto-generated destructor stub
}

