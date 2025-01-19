/*
 * Rays.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Rays.h"

Rays::Rays(int nRays, std::vector<double> sourcePosition) : numRays(nRays), sourcePosition(sourcePosition) {
	direction = std::vector<std::vector<double>>(nRays, std::vector<double>(3, 0.0));
	initializeDirections();
	position = std::vector<std::vector<double>>(nRays, sourcePosition);

	opticalDepth = std::vector<double>(nRays, 0.0);
	numTraversedCells = std::vector<double>(nRays, 0.0);
	insideDomain = std::vector<bool>(nRays, true);

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

/*
void Rays::doRayTracing(){
	for(int iRay = 0; iRay < numRays; iRay++){

		int iCell = startCell;
		while(insideDomain[iRay]){

			if(boundaryFlag[iCell]){
				exitCell, distanceToExit = getDistanceToBoundary(domain, rays, iRay);
				insideDomain[iRay] = false;
			}
			else{
				exitCell, distanceToExit = findNextCell(iCell);
			}

		}

	}

}
*/

Rays::~Rays() {
	// TODO Auto-generated destructor stub
}

