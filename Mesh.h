/*
 * Mesh.h
 *
 *  Created on: 19 Jan 2025
 *      Author: ntc132
 */

#ifndef MESH_H_
#define MESH_H_

class Mesh {
public:
	Mesh(std::string fileMeshIndices);
	virtual ~Mesh();

	int numCells;
	std::vector<int> cellIDs;
	std::vector<std::vector<int>> neighbourList;

	std::vector<std::vector<int> > collectNeighbours(const std::vector<std::pair<int, int>>& IdPairs, std::vector<int> cellIDs);
	std::vector<int> getCellIDs(const std::vector<std::pair<int, int>>& IdPairs);
	std::vector<std::pair<int, int>> readVoronoiIndices(const std::string& filename);

protected:
	std::string fileMeshIndices;
	std::vector<std::pair<int, int>> IdPairs;

};

#endif /* MESH_H_ */
