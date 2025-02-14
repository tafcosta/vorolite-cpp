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
	startCell = mesh.findHostCellID(sourcePosition, -1)[0];

	std::vector<float> floatCoordinates = mesh.cellCoordinates[startCell];
	std::vector<double> doubleCoordinates;
	doubleCoordinates.reserve(floatCoordinates.size());

	for (float value : floatCoordinates)
	    doubleCoordinates.push_back(static_cast<double>(value));





	rayDirection = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	theta = std::vector<double>(nRays, 0.0);
	phi = std::vector<double>(nRays, 0.0);

	visitedCells = std::vector<std::vector<int>>(nRays);

	rayPosition = std::vector<std::vector<double>>(nRays, doubleCoordinates);
	initializeDirections();

	columnDensity     = std::vector<double>(nRays, 0.0);
	distanceTravelled = std::vector<double>(nRays, 0.0);
	insideDomain      = std::vector<bool>(nRays, true);
	ignoreRay         = std::vector<bool>(nRays, false);
	numTraversedCells = std::vector<double>(nRays, 0.0);

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

 int Rays::findNextCell(int iCell, int iRay, bool verbose){

	 double distanceToExit = std::numeric_limits<double>::max();
	 double distanceToExitTmp;
	 double overshoot; // new

	 std::ostringstream debugOutput;
	 int exitCell  = -1;

	 std::vector<int> closestCells;
	 std::vector<float> cellPos = mesh.cellCoordinates[iCell];
	 std::vector<double> normalVector (3, 0.0);
	 std::vector<double> pointOnInterface (3, 0.0);
	 std::vector<double> positionTmp (3, 0.0); // new

	 if(verbose){
		 double distanceBetRayAndCell = sqrt((cellPos[0] - rayPosition[iRay][0]) * (cellPos[0] - rayPosition[iRay][0]) + (cellPos[1] - rayPosition[iRay][1]) * (cellPos[1] - rayPosition[iRay][1]) + (cellPos[2] - rayPosition[iRay][2]) * (cellPos[2] - rayPosition[iRay][2]));
		 debugOutput << "Host Index = " << iCell << " Host Cell Position = " << cellPos[0] << ", " << cellPos[1] << ", " << cellPos[2] << " Distance from Ray to Cell (before update) = " << distanceBetRayAndCell << "\n";
	 }


	 for (int neighbour : mesh.neighbourList[iCell]){
		 std::vector<float> neighbourPos = mesh.cellCoordinates[neighbour];

		 if(verbose){
			 debugOutput << "Neighbour Index = " << neighbour << ", Neighbour Position = " << neighbourPos[0] << ", " << neighbourPos[1] << ", " << neighbourPos[2] << "\n";
		 }


		 for (int i = 0; i < 3; i++){
			 normalVector[i] = neighbourPos[i] - cellPos[i];
			 pointOnInterface[i] = 0.5 * (cellPos[i] + neighbourPos[i]);
		 }

		 double denominator = normalVector[0] * rayDirection[iRay][0] + normalVector[1] * rayDirection[iRay][1] + normalVector[2] * rayDirection[iRay][2];

		 distanceToExitTmp = normalVector[0] * (pointOnInterface[0] - rayPosition[iRay][0])
				 + normalVector[1] * (pointOnInterface[1] - rayPosition[iRay][1])
				 + normalVector[2] * (pointOnInterface[2] - rayPosition[iRay][2]);

		 if(denominator <= 0.)
			 continue;

		 distanceToExitTmp = distanceToExitTmp/denominator;

		 if(verbose)
			 debugOutput << std::scientific << std::setprecision(12) << "Distance to Neighbour = " << distanceToExitTmp << " denominator = " << denominator << "\n";


		 if((distanceToExitTmp <= distanceToExit) && (distanceToExitTmp > 0)){
			 // there was another test here that I have removed as it should not be needed anymore; however I was not sure whether it was correct in the original version of the function

			 distanceToExit = distanceToExitTmp;
			 exitCell = neighbour;

			 if(verbose)
				 debugOutput << "New Neighbour candidate = " << neighbour << "\n";
		 }
	 }


	 if(distanceToExit > mesh.boxSize){
	 	insideDomain[iRay] = false;
	 	distanceToExit = 0.0;
	 	exitCell = -1;
	 	std::cout << "distance to exit larger than box size; stopped!" << "\n";
	 }

	 if (distanceToExit < 1e-10)
	 	std::cout << "distanceToExit = " << distanceToExit << "You seem to have ended up on an edge, how did you do that?!" << "\n";

	columnDensity[iRay] += distanceToExit * mesh.cellDensity[iCell]; // we add up whatever density we have encountered on the way out of the cell

	// we now know where we would pierce the face of the cell; we want to continue onwards

	overshoot = distanceToExit / 10;
	do {
		for (int i = 0; i < 3; i++)
			 positionTmp[i] = rayPosition[iRay][i] + rayDirection[iRay][i] * (distanceToExit + overshoot);

		overshoot /= 2;

	} while (mesh.findHostCellID(positionTmp, -1)[0] != exitCell);  // todo maybe replace by a check of the entire list, but a priori there should only be one cell as we no longer are on an edge





	 if(verbose){
		 debugOutput << "Ray position (before update) = " << rayPosition[iRay][0] << ", "
				 << rayPosition[iRay][1] << ", " << rayPosition[iRay][2] << "\n";
	    }


	 for (int i = 0; i < 3; i++)
		 rayPosition[iRay][i] += rayDirection[iRay][i] * (distanceToExit + overshoot);


	 visitedCells[iRay].push_back(iCell);
	 columnDensity[iRay] += overshoot * mesh.cellDensity[exitCell]; // we add up all the density we encounter in the next cell up to the new position
	 distanceTravelled[iRay] += distanceToExit + overshoot;
	 numTraversedCells[iRay] += 1;

	 closestCells = mesh.findHostCellID(positionTmp, -1);

	 if(verbose){
		 double distanceToCell = mesh.getDistanceToCell(rayPosition[iRay], exitCell);

		 debugOutput << "Ray position (after update) = " << rayPosition[iRay][0] << " "
				 << rayPosition[iRay][1] << " "
				 << rayPosition[iRay][2] << "\n";

		 debugOutput  << "Closest cells = ";

		 for(int i = 0; i < closestCells.size(); i++)
			 debugOutput << closestCells[i] << ", ";

		 debugOutput  << "\n";

		 debugOutput  << "Distance from ray to exit cell = " << distanceToCell << "\n";

		 if(exitCell != -1){
			 std::vector<float> cellPos = mesh.cellCoordinates[exitCell];
			 debugOutput << "Next cell (from neighbour) = " << exitCell << "\n"
					 << "Position of cell (from neighbour search) = "
					 << cellPos[0] << ", " << cellPos[1] << ", " << cellPos[2] << "\n";
		 }
	 }


	 if(rayPosition[iRay][0] > mesh.boxSize || rayPosition[iRay][0] < 0 || rayPosition[iRay][1] > mesh.boxSize || rayPosition[iRay][1] < 0 || rayPosition[iRay][2] > mesh.boxSize || rayPosition[iRay][2] < 0 || distanceTravelled[iRay] >= maxRadius){
		 insideDomain[iRay] = false;
		 exitCell = -1;
		std::cout << "new position outside box; stopped" << "\n";
	 } else {

		 bool contains = std::find(closestCells.begin(), closestCells.end(), exitCell) != closestCells.end();

		 if(exitCell == -1){ //contains == false ||
			 warningIssued = true;

			 /*
			 if(verbose == false)
				 std::cerr << "Warning: either inside domain and no neighbours found or mismatch in ray and neigbour cells. Ray number = " << iRay << ", exitCell = " << exitCell << std::endl;
			*/

			 if(verbose)
				 std::cout << debugOutput.str() << std::endl;  // Output all debug information at once

			 ignoreRay[iRay] = true;
		 }
	 }

	 return exitCell;
 }



	 void Rays::outputResults() {
     std::ofstream outputFile("ray_output.txt");

     if (!outputFile.is_open()) {
         std::cerr << "Error opening file for output!" << std::endl;
         return;
     }

     outputFile << "Ray\tTheta\tPhi\tColumn\tIgnore" << std::endl;

     for (int i = 0; i < numRays; ++i) {
         outputFile << i << "\t"  // Ray number
                    << std::fixed << std::setprecision(3) << theta[i] << "\t"
                    << std::fixed << std::setprecision(3) << phi[i] << "\t"
                    << std::fixed << std::setprecision(6) << columnDensity[i] << "\t"
					<< std::fixed << ignoreRay[i] << std::endl;

     }

     outputFile.close();
     std::cout << "Results have been written to 'ray_output.txt'" << std::endl;


     /*
     std::ofstream visitedCellsFile("ray_visited_cells.txt");

     if (!visitedCellsFile.is_open()) {
         std::cerr << "Error opening file for visited cells output!" << std::endl;
         return;
     }

     visitedCellsFile << "Ray\tVisitedCells" << std::endl;

     for (int i = 0; i < numRays; ++i) {
         visitedCellsFile << i << "\t";  // Ray number

         // Output the visited cells for this ray
         for (int visitedCell : visitedCells[i]) {
             visitedCellsFile << visitedCell << " ";
         }
         visitedCellsFile << std::endl;
     }

     visitedCellsFile.close();
     */

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
			iCell = findNextCell(iCellOld, iRay, false);

			if(debug){
				if(warningIssued){

					for(int i = 0; i < 3; i++)
						rayPosition[iRay][i] = oldRayPosition[i];

					iCell = findNextCell(iCellOld, iRay, true);
					insideDomain[iRay] = false;
				}
			}
		}
	}
}


Rays::~Rays() {
	// TODO Auto-generated destructor stub
}

