/*
 * Mesh.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Mesh.h"

Mesh::Mesh(std::string fileMeshIndices, std::string snapshot) : fileMeshIndices(fileMeshIndices), snapshot(snapshot) {

    readSnapshot(snapshot);
	isAtBoundary.resize(numCells, false);

	IdPairs       = readVoronoiIndices(fileMeshIndices);
	neighbourList = collectNeighbours(IdPairs, cellIDs);
}

void Mesh::readSnapshot(const std::string& snapshot) {
    try {
        H5::H5File file(snapshot, H5F_ACC_RDONLY);

        H5::DataSet densityDataset = file.openDataSet("/PartType0/Density");
        H5::DataSet coordinatesDataset = file.openDataSet("/PartType0/Coordinates");
        H5::DataSet idDataset = file.openDataSet("/PartType0/ParticleIDs");

        H5::Group headerGroup = file.openGroup("/Header");
        H5::Attribute boxSizeAttribute = headerGroup.openAttribute("BoxSize");
        boxSizeAttribute.read(H5::PredType::NATIVE_DOUBLE, &boxSize);

        H5::DataSpace densitySpace = densityDataset.getSpace();
        hsize_t numDensities;
        densitySpace.getSimpleExtentDims(&numDensities);
        cellDensity.resize(numDensities);
        densityDataset.read(cellDensity.data(), H5::PredType::NATIVE_DOUBLE);
        numCells = cellDensity.size();

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

        H5::DataSpace idSpace = idDataset.getSpace();
        hsize_t numIDs;
        idSpace.getSimpleExtentDims(&numIDs);
        cellIDs.resize(numIDs);
        idDataset.read(cellIDs.data(), H5::PredType::NATIVE_INT);

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


int Mesh::findHostCellID(const std::vector<double>& target, int cellSkip) {
    int closestCell;
    double minDistance = std::numeric_limits<double>::infinity();

    for (int iCell = 0; iCell < numCells; iCell++) {

    	if(iCell == cellSkip)
    		continue;

        double dist = squaredDistance(cellCoordinates[iCell], target);

        if (dist < minDistance) {
            minDistance = dist;
            closestCell = iCell;
        }
    }

    return closestCell;
}


Mesh::~Mesh() {

}

