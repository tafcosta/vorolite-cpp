/*
 * Mesh.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Mesh.h"

Mesh::Mesh(std::string fileMeshIndices, std::string snapshot) : fileMeshIndices(fileMeshIndices), snapshot(snapshot) {

	IdPairs = readVoronoiIndices(fileMeshIndices);
	cellIDs = getCellIDs(IdPairs);
	numCells = cellIDs.size();

	isAtBoundary.resize(numCells, false);

	neighbourList = collectNeighbours(IdPairs, cellIDs);

    readSnapshot(snapshot);
}

void Mesh::readSnapshot(const std::string& snapshot) {
    try {
        std::cout << "Reading Arepo snapshot..." << std::endl;

        H5::H5File file(snapshot, H5F_ACC_RDONLY);

        H5::DataSet densityDataset = file.openDataSet("/PartType0/Density");
        H5::DataSet coordinatesDataset = file.openDataSet("/PartType0/Coordinates");

        H5::DataSpace densitySpace = densityDataset.getSpace();
        hsize_t numDensities;
        densitySpace.getSimpleExtentDims(&numDensities);
        cellDensity.resize(numDensities);
        densityDataset.read(cellDensity.data(), H5::PredType::NATIVE_DOUBLE);

        if (numCells != cellDensity.size()) {
            std::cerr << "Error: numCells (" << numCells << ") does not match the number of cells from mesh output... (" << cellDensity.size() << ")" << std::endl;
            exit(EXIT_FAILURE);
        }

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

        std::cout << "Done reading snapshot. There are " << cellDensity.size() << " cells." << std::endl;

    } catch (H5::Exception& e) {
        std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
    }
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

std::vector<int> Mesh::getCellIDs(const std::vector<std::pair<int, int>>& IdPairs) {
  std::unordered_set<int> uniqueIDs;

  for (const auto& pair : IdPairs)
    uniqueIDs.insert(pair.first);

  return std::vector<int> (uniqueIDs.begin(), uniqueIDs.end());
}

std::vector<std::vector<int> > Mesh::collectNeighbours(const std::vector<std::pair<int, int>>& IdPairs, std::vector<int> cellIDs) {
    std::unordered_map<int, std::vector<int>> neighboursMap;
    std::vector<std::vector<int>> neighbourList;

    for (const auto& pair : IdPairs)
      neighboursMap[pair.first].push_back(pair.second);

    for (const int id : cellIDs) {
        std::vector<int> row;

        if (neighboursMap.find(id) != neighboursMap.end())
            row = neighboursMap[id];

        neighbourList.push_back(row);
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


int Mesh::findHostCellID(const std::vector<double>& target) {
    int closestCellID;
    double minDistance = std::numeric_limits<double>::infinity();


    for (int iCell = 0; iCell < numCells; iCell++) {

        double dist = squaredDistance(cellCoordinates[iCell], target);

        /*
        if (dist < minDistance) {
            minDistance = dist;
            closestCellID = cellIDs[iCell];
        }
        */

    }

    return closestCellID;
}


Mesh::~Mesh() {

}

