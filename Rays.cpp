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
	direction = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	initializeDirections();
	position = std::vector<std::vector<double>>(nRays, sourcePosition);

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

        direction[iRay][0] = std::cos(phi) * std::sin(theta);
        direction[iRay][1] = std::sin(phi) * std::sin(theta);
        direction[iRay][2] = std::cos(theta);
    }
}

 void Rays::findNextCell(int iCell){

	 double distanceToExit = std::numeric_limits<double>::max();
	 int exitCell  = -1;

	 std::vector<float> cellPos = mesh.cellCoordinates[iCell];

	 for (int neighbour : mesh.neighbourList[iCell]){
		 std::vector<float> neighbourPos = mesh.cellCoordinates[neighbour];
		 std::vector<float> normalVector (3, 0.0);

	        if (cellPos.size() == neighbourPos.size()) {

	        	normalVector[0] = cellPos[0] - neighbourPos[0];
	        	normalVector[1] = cellPos[1] - neighbourPos[1];
	        	normalVector[2] = cellPos[2] - neighbourPos[2];

	            double distance = std::sqrt(normalVector[0] * normalVector[0] + normalVector[1] * normalVector[1] + normalVector[2] * normalVector[2]);
	            std::cout << "cell pos: " << cellPos[1] << ", neighbour pos: " << neighbourPos[1] << std::endl;
	        }
	        else {
	            std::cerr << "Error: Inconsistent vector sizes!" << std::endl;
	        }

	 }

	 /*
    for iNeighbour in range(numberNeighbours[iCell]):
        normalVector       = pointCloud[iCell] - pointCloud[indices[indptr[iCell]:indptr[iCell + 1]][iNeighbour]]
        pointOnInterface   = (pointCloud[iCell] + pointCloud[indices[indptr[iCell]:indptr[iCell + 1]][iNeighbour]]) / 2.0
        denominator        = np.dot(normalVector, np.array([rays.kx[iRay], rays.ky[iRay], rays.kz[iRay]]).ravel())

        if(np.abs(denominator) > sys.float_info.epsilon):
            distanceToExitTmp  = np.dot(normalVector,  (pointOnInterface - np.array([rays.xPos[iRay], rays.yPos[iRay], rays.zPos[iRay]]).ravel())) / denominator
        else:
            continue

        if ((distanceToExitTmp < distanceToExit) and (distanceToExitTmp > 0)):
            distanceToExit = distanceToExitTmp
            exitCell = indices[indptr[iCell]:indptr[iCell + 1]][iNeighbour]
            numberPossibleNeighbours += 1

    if(numberPossibleNeighbours == 0):
        raise ValueError(f"Error: Loop terminated because no elligible neighbours were found.")

    if(exitCell == -1):
        raise ValueError(f"No exit cell found.")

    */
 }


void Rays::doRayTracing(){


	for(int iRay = 0; iRay < numRays; iRay++){

		int iCell = startCell;
		findNextCell(iCell);

		/*
		while(insideDomain[iRay]){

			if(mesh.isAtBoundary[iCell]){
				//getDistanceToBoundary(domain, rays, iRay);
				insideDomain[iRay] = false;
			}
			else{
				findNextCell(iCell);
				insideDomain[iRay] = false;
			}
		}
		*/

	}


	std::cout << "Do nothing" << std::endl;


}


Rays::~Rays() {
	// TODO Auto-generated destructor stub
}

