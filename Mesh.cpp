/*
 * Mesh.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Mesh.h"
#include <filesystem>

Mesh::Mesh(std::string fileMeshIndices, std::string snapshot, double maxRadius, std::vector<double> sourcePosition) : fileMeshIndices(fileMeshIndices), snapshot(snapshot), maxRadius(maxRadius), sourcePosition(sourcePosition) {

    readSnapshot(snapshot);
	getNumCellsInRegion();

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

double Mesh::getIncomingFlux(int iCell){
	return cellIncomingFlux[iCell];
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


void Mesh::readSnapshot(const std::string& snapshotBase) {
    try {
    	std::vector<std::string> files = getSnapshotFiles(snapshotBase);

    	std::cout << "[DEBUG] getSnapshotFiles() returned " << files.size() << " files:" << std::endl;
    	for (const auto& f : files) {
    	    std::cout << " - " << f << std::endl;
    	}

        bool headerRead = false;

        for (const std::string& fileName : files) {
            H5::H5File file(fileName, H5F_ACC_RDONLY);

            if (!headerRead) {
                readHeader(file);
                headerRead = true;
            }

            appendDensity(file);
            appendMass(file);
            appendIDs(file);
            appendCoordinates(file);
            appendVelocities(file);
        }

        numCells = cellDensity.size();
        cellIndices.resize(numCells);
        std::iota(cellIndices.begin(), cellIndices.end(), 0);

        std::cout << "Done reading snapshot: " << snapshotBase << std::endl;
        std::cout << "Total number of cells: " << numCells << std::endl;

    } catch (H5::Exception& e) {
        std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
    }

    // TODO: cache IC data if needed
}

void Mesh::readHeader(H5::H5File& file) {
    H5::Group headerGroup = file.openGroup("/Header");

    headerGroup.openAttribute("BoxSize").read(H5::PredType::NATIVE_DOUBLE, &boxSize);
    headerGroup.openAttribute("UnitLength_in_cm").read(H5::PredType::NATIVE_DOUBLE, &unitLength);
    headerGroup.openAttribute("UnitMass_in_g").read(H5::PredType::NATIVE_DOUBLE, &unitMass);
    headerGroup.openAttribute("UnitVelocity_in_cm_per_s").read(H5::PredType::NATIVE_DOUBLE, &unitVelocity);
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

void Mesh::appendDensity(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/Density");
    H5::DataSpace space = dataset.getSpace();

    hsize_t numElements;
    space.getSimpleExtentDims(&numElements);

    std::vector<double> buffer(numElements);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    cellDensity.insert(cellDensity.end(), buffer.begin(), buffer.end());
}

void Mesh::appendMass(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/Masses");
    H5::DataSpace space = dataset.getSpace();

    hsize_t numElements;
    space.getSimpleExtentDims(&numElements);

    std::vector<double> buffer(numElements);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    cellMass.insert(cellMass.end(), buffer.begin(), buffer.end());
}

void Mesh::appendIDs(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/ParticleIDs");
    H5::DataSpace space = dataset.getSpace();

    hsize_t numElements;
    space.getSimpleExtentDims(&numElements);

    std::vector<int> buffer(numElements);
    dataset.read(buffer.data(), H5::PredType::NATIVE_INT);
    cellIDs.insert(cellIDs.end(), buffer.begin(), buffer.end());
}


void Mesh::appendCoordinates(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/Coordinates");
    H5::DataSpace space = dataset.getSpace();

    hsize_t dims[2];
    space.getSimpleExtentDims(dims);

    std::vector<float> buffer(dims[0] * dims[1]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_FLOAT);

    for (hsize_t i = 0; i < dims[0]; ++i) {
        std::vector<float> row(dims[1]);
        for (hsize_t j = 0; j < dims[1]; ++j) {
            row[j] = buffer[i * dims[1] + j];
        }
        cellCoordinates.push_back(std::move(row));
    }
}

void Mesh::appendVelocities(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/Velocities");
    H5::DataSpace space = dataset.getSpace();

    hsize_t dims[2];
    space.getSimpleExtentDims(dims);

    std::vector<float> buffer(dims[0] * dims[1]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_FLOAT);

    for (hsize_t i = 0; i < dims[0]; ++i) {
        std::vector<float> row(dims[1]);
        for (hsize_t j = 0; j < dims[1]; ++j) {
            row[j] = buffer[i * dims[1] + j];
        }
        cellVelocities.push_back(std::move(row));
    }
}


std::vector<std::string> Mesh::getSnapshotFiles(const std::string& snapshotPath) {
    std::vector<std::string> files;

    std::filesystem::path inputPath(snapshotPath);
    std::filesystem::path dir = inputPath.parent_path();
    if (dir.empty()) dir = ".";

    std::string baseName = inputPath.stem().string();     // "snap_006"
    std::string extension = inputPath.extension().string(); // ".hdf5"

    // Step 1: Look for split files like snap_006.0.hdf5, snap_006.1.hdf5, ...
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;

        std::string fname = entry.path().filename().string();

        // Match baseName.N.hdf5, e.g., snap_006.0.hdf5
        if (fname.rfind(baseName + ".", 0) == 0 &&
            fname.size() > baseName.size() + 6 &&
            fname.compare(fname.size() - 6, 6, ".hdf5") == 0) {
            files.push_back(entry.path().string());
        }
    }

    // Step 2: Fallback to monolithic file if no split files found
    if (files.empty() && std::filesystem::exists(snapshotPath)) {
        files.push_back(snapshotPath);
    }

    std::sort(files.begin(), files.end());
    return files;
}



Mesh::~Mesh() {

}

