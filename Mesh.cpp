/*
 * Mesh.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Mesh.h"

Mesh::Mesh(std::string fileMeshIndices, std::string snapshot, double maxRadius, std::vector<double> sourcePosition) : fileMeshIndices(fileMeshIndices), snapshot(snapshot), maxRadius(maxRadius), sourcePosition(sourcePosition) {

    readSnapshot(snapshot);
	getNumCellsInRegion();

	for(int iCell = 0; iCell < numCells; iCell++){
		cellDensity[iCell] = cellDensity[iCell]/10;
		cellMass[iCell] = cellMass[iCell]/10;
	}


	cellVisitsByRay.resize(numCells, 0);
    cellFlux.resize(numCells, 0.0);
    cellIncomingFlux.resize(numCells, 0.0);
    cellLocalColumn.resize(numCells, 0.0);
	cellHIIFraction.resize(numCells, 0.0);
	fluxOfRayInCell.resize(numCells); //The first dimension should be number of rays

	IdPairs       = readVoronoiIndices(fileMeshIndices);
	neighbourList = collectNeighbours(IdPairs, cellIDs);
}

void Mesh::resizeFluxOfRayInCell(int iRay, int numVisitedCells){
	fluxOfRayInCell[iRay].resize(numVisitedCells, 0.0);
}

double Mesh::getMass(int iCell){
	return cellMass[iCell];
}

double Mesh::getDensity(int iCell){
	return cellDensity[iCell];
}

double Mesh::getNumberDensity_in_cgs(int iCell){
	return cellDensity[iCell] * unitMass / unitLength / unitLength / unitLength / protonMass;
}

double Mesh::getElectronNumberDensity_in_cgs(int iCell){
	return cellHIIFraction[iCell] * getNumberDensity_in_cgs(iCell);
}

double Mesh::getMeanMolecularWeight(int iCell){
	return 1./(1 + cellHIIFraction[iCell]);
}

double Mesh::getFluxOfRayInCell(int iRay, int iCell){
	return fluxOfRayInCell[iRay][iCell];
}

void Mesh::setFluxOfRayInCell(int iRay, int iCell, double newValue){
	fluxOfRayInCell[iRay][iCell] = newValue;
}

void Mesh::setHIIFraction(int iCell, double newValue){
	cellHIIFraction[iCell] = newValue;

	if(newValue > 1)
		cellHIIFraction[iCell] = 1;

	if(newValue < 1e-5)
		cellHIIFraction[iCell] = 1.e-5;
}

double Mesh::getFlux(int iCell){
	return cellFlux[iCell];
}

double Mesh::getHIIFraction(int iCell){
	return cellHIIFraction[iCell];
}

int Mesh::getIndex(int iCell){
	return cellIndices[iCell];
}

void Mesh::getNumCellsInRegion(){

    std::vector<std::vector<float>> filteredCoordinates;
    std::vector<std::vector<float>> filteredVelocities;
    std::vector<double> filteredDensity;
    std::vector<double> filteredMasses;
    std::vector<int> filteredIDs;
    std::vector<int> filteredCellIndices;

	std::vector<float> cellPos;

	int startCell = findHostCellID(sourcePosition, -1)[0];

	for(int iCell = 0; iCell < numCells; iCell++){
		cellPos   = cellCoordinates[iCell];
        double dx = cellPos[0] - cellCoordinates[startCell][0]; //sourcePosition[0];
        double dy = cellPos[1] - cellCoordinates[startCell][1]; //sourcePosition[1];
        double dz = cellPos[2] - cellCoordinates[startCell][2]; //sourcePosition[2];
        double rDistance = std::sqrt(dx*dx + dy*dy + dz*dz);

		if(rDistance <= 1.2 * maxRadius){
            filteredCoordinates.push_back(cellCoordinates[iCell]);
            filteredVelocities.push_back(cellVelocities[iCell]);
            filteredDensity.push_back(cellDensity[iCell]);
            filteredMasses.push_back(cellMass[iCell]);
            filteredIDs.push_back(cellIDs[iCell]);
            filteredCellIndices.push_back(cellIndices[iCell]);
		}
	}

    cellCoordinates = std::move(filteredCoordinates);
    cellVelocities = std::move(filteredVelocities);
    cellDensity = std::move(filteredDensity);
    cellIDs = std::move(filteredIDs);
    cellIndices = std::move(filteredCellIndices);
    cellMass = std::move(filteredMasses);

    numCells = cellDensity.size();

    std::cout << "Reduced to " << numCells << " cells within maxRadius = " << maxRadius << std::endl;
}

void Mesh::resetFluxes(){
	for(int iCell = 0; iCell < numCells; iCell++){
		cellFlux[iCell] = 0.;
		cellIncomingFlux[iCell] = 0.;
	}
}

void Mesh::readSnapshot(const std::string& snapshot) {
    try {
        H5::H5File file(snapshot, H5F_ACC_RDONLY);

        H5::DataSet densityDataset = file.openDataSet("/PartType0/Density");
        H5::DataSet massDataset = file.openDataSet("/PartType0/Masses");
        H5::DataSet coordinatesDataset = file.openDataSet("/PartType0/Coordinates");
        H5::DataSet velocitiesDataset = file.openDataSet("/PartType0/Velocities");
        H5::DataSet idDataset = file.openDataSet("/PartType0/ParticleIDs");

        H5::Group headerGroup = file.openGroup("/Header");
        H5::Attribute boxSizeAttribute = headerGroup.openAttribute("BoxSize");
        boxSizeAttribute.read(H5::PredType::NATIVE_DOUBLE, &boxSize);

        H5::Attribute attrLength = headerGroup.openAttribute("UnitLength_in_cm");
        attrLength.read(H5::PredType::NATIVE_DOUBLE, &unitLength);
        H5::Attribute attrMass = headerGroup.openAttribute("UnitMass_in_g");
        attrMass.read(H5::PredType::NATIVE_DOUBLE, &unitMass);
        H5::Attribute attrVelocity = headerGroup.openAttribute("UnitVelocity_in_cm_per_s");
        attrVelocity.read(H5::PredType::NATIVE_DOUBLE, &unitVelocity);

        H5::DataSpace densitySpace = densityDataset.getSpace();
        hsize_t numDensities;
        densitySpace.getSimpleExtentDims(&numDensities);
        cellDensity.resize(numDensities);
        densityDataset.read(cellDensity.data(), H5::PredType::NATIVE_DOUBLE);
        numCells = cellDensity.size();

        H5::DataSpace massSpace = massDataset.getSpace();
        hsize_t numMass;
        massSpace.getSimpleExtentDims(&numMass);
        cellMass.resize(numMass);
        massDataset.read(cellMass.data(), H5::PredType::NATIVE_DOUBLE);

    	cellIndices.resize(numDensities, 0);
    	std::iota(cellIndices.begin(), cellIndices.end(), 0);

        H5::DataSpace coordinatesSpace = coordinatesDataset.getSpace();
        hsize_t dims[2];
        coordinatesSpace.getSimpleExtentDims(dims);
        std::vector<float> cellCoordinates1D(dims[0] * dims[1]);
        coordinatesDataset.read(cellCoordinates1D.data(), H5::PredType::NATIVE_FLOAT);

        cellCoordinates.resize(dims[0], std::vector<float>(dims[1]));

        for (size_t row = 0; row < dims[0]; ++row) {
             for (size_t col = 0; col < dims[1]; ++col) {
                 size_t index = row * dims[1] + col;
                 cellCoordinates[row][col] = cellCoordinates1D[index];
             }
         }


        H5::DataSpace velocitiesSpace = velocitiesDataset.getSpace();
        velocitiesSpace.getSimpleExtentDims(dims);
        std::vector<float> cellVelocities1D(dims[0] * dims[1]);
        velocitiesDataset.read(cellVelocities1D.data(), H5::PredType::NATIVE_FLOAT);

        cellVelocities.resize(dims[0], std::vector<float>(dims[1]));

        for (size_t row = 0; row < dims[0]; ++row) {
             for (size_t col = 0; col < dims[1]; ++col) {
                 size_t index = row * dims[1] + col;
                 cellVelocities[row][col] = cellVelocities1D[index];
             }
         }

        H5::DataSpace idSpace = idDataset.getSpace();
        hsize_t numIDs;
        idSpace.getSimpleExtentDims(&numIDs);
        cellIDs.resize(numIDs);
        idDataset.read(cellIDs.data(), H5::PredType::NATIVE_INT);

        std::cout << "Done reading snapshot." << std::endl;
        std::cout << "There are " << cellDensity.size() << " cells." << std::endl;

    } catch (H5::Exception& e) {
        std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
    }

    //todo: save relevant IC data to prevent reading from snapshot every time.
}

std::vector<std::pair<int, int>> Mesh::readVoronoiIndices(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return {};
    }

    // Get file size
    file.seekg(0, std::ios::end);
    size_t fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    // Read entire file into buffer
    std::vector<char> buffer(fileSize);
    file.read(buffer.data(), fileSize);
    file.close();

    const char* data = buffer.data();
    const char* end = data + fileSize;

    // Skip initial offset
    data += sizeof(int);

    std::vector<std::pair<int, int>> IdPairs;

    while (data + 3 * sizeof(int) <= end) {
        int faceId      = *reinterpret_cast<const int*>(data);
        int faceIdOther = *reinterpret_cast<const int*>(data + sizeof(int));
        int numVertices = *reinterpret_cast<const int*>(data + 2 * sizeof(int));
        data += 3 * sizeof(int);

        IdPairs.emplace_back(faceId, faceIdOther);

        // Skip vertices
        data += numVertices * sizeof(int);
    }

    std::cout << "Collected " << IdPairs.size() << " cell pairs..." << std::endl;
    return IdPairs;
}

std::vector<std::vector<int>> Mesh::collectNeighbours(const std::vector<std::pair<int, int>>& IdPairs, std::vector<int>& cellIDs)
{
    std::unordered_map<int, int> cellIDToIndex;
    std::unordered_map<int, std::unordered_set<int>> neighboursMap;
    std::vector<std::vector<int>> neighbourList;

    for (size_t i = 0; i < cellIDs.size(); ++i)
        cellIDToIndex[cellIDs[i]] = i;

    for (const auto& pair : IdPairs) {
        int cell1 = pair.first;
        int cell2 = pair.second;

        if (cellIDToIndex.find(cell1) != cellIDToIndex.end() && cellIDToIndex.find(cell2) != cellIDToIndex.end()) {
            int index1 = cellIDToIndex[cell1];
            int index2 = cellIDToIndex[cell2];

            neighboursMap[index1].insert(index2);
            neighboursMap[index2].insert(index1);
        }
    }

    for (size_t i = 0; i < cellIDs.size(); ++i) {
        if (neighboursMap.find(i) != neighboursMap.end()) {
            neighbourList.push_back(std::vector<int>(neighboursMap[i].begin(), neighboursMap[i].end()));
        } else {
            neighbourList.push_back({});
        }
    }

    return neighbourList;
}


double Mesh::squaredDistance(const std::vector<float>& point1, const std::vector<double>& point2) {
    double dist = 0.0;
    for (size_t i = 0; i < point1.size(); ++i) {
        dist += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    return dist;
}

std::vector<int> Mesh::findHostCellID(const std::vector<double>& target, int cellGuess) {
    std::vector<int> closestCells;
    std::vector<int> possibleCells;
    closestCells.reserve(numCells);
    possibleCells.reserve(numCells);

    if (cellGuess != -1) {
    	possibleCells.push_back(cellGuess);
    	for (int neighbour : neighbourList[cellGuess])
    		possibleCells.push_back(neighbour);
    }
    else{
    	possibleCells.resize(numCells);
        std::iota(possibleCells.begin(), possibleCells.end(), 0);
    }

    double minDistance = std::numeric_limits<double>::infinity();

    for (int iCell : possibleCells) {
        double distance = squaredDistance(cellCoordinates[iCell], target);

        if (std::fabs(distance - minDistance) < 1.e-10) {
            closestCells.push_back(iCell);
        } else if (distance < minDistance) {
            closestCells = {iCell};
            minDistance = distance;
        }
    }

    return closestCells;
}


bool Mesh::checkIfExitCellNeighboursCurrentCell(int iCell, int exitCell){
    bool test = false;
    std::vector<int> possibleCells;

	possibleCells.push_back(iCell);
    possibleCells.insert(possibleCells.end(), neighbourList[iCell].begin(), neighbourList[iCell].end());

    for (int cell : possibleCells) {
    	if(cell == exitCell)
    		test = true;
    }

    return test;
}

double Mesh::getDistanceToCell(const std::vector<double>& target, int cellIndex) {
    return sqrt(squaredDistance(cellCoordinates[cellIndex], target));

}

double Mesh::getDistanceBetweenCells(int iCell, int jCell) {

	double dist = (cellCoordinates[iCell][0] - cellCoordinates[jCell][0]) * (cellCoordinates[iCell][0] - cellCoordinates[jCell][0])
        		+ (cellCoordinates[iCell][1] - cellCoordinates[jCell][1]) * (cellCoordinates[iCell][1] - cellCoordinates[jCell][1])
        		+ (cellCoordinates[iCell][2] - cellCoordinates[jCell][2]) * (cellCoordinates[iCell][2] - cellCoordinates[jCell][2]);

	dist = sqrt(dist);

    return dist;

}

void Mesh::saveVoronoiIndices(const std::string& filename, const std::vector<std::pair<int, int>>& IdPairs) {
    std::ofstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    size_t size = IdPairs.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));

    for (const auto& pair : IdPairs) {
        file.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
        file.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second));
    }

    file.close();
    std::cout << "Saved " << IdPairs.size() << " cell pairs to " << filename << std::endl;
}


Mesh::~Mesh() {

}

