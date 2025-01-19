/*
 * Mesh.cpp
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#include "common_includes.h"
#include "Mesh.h"

Mesh::Mesh(std::string fileMeshIndices) : fileMeshIndices(fileMeshIndices) {

	IdPairs = readVoronoiIndices(fileMeshIndices);
	cellIDs = getCellIDs(IdPairs);
	neighbourList = collectNeighbours(IdPairs, cellIDs);

	numCells = cellIDs.size();

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


Mesh::~Mesh() {

}

