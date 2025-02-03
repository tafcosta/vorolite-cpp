/*
 * Mesh.h
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#ifndef MESH_H_
#define MESH_H_

class Mesh {
public:
	Mesh(std::string fileMeshIndices, std::string snapshot);
	virtual ~Mesh();

	int numCells;
	std::vector<int> cellIDs;
    std::vector<std::vector<float>> cellCoordinates;
    std::vector<double> cellDensity;
    std::vector<bool> isAtBoundary;

	std::vector<std::vector<int>> neighbourList;
	int findHostCellID(const std::vector<double>& target);


private:
	std::vector<std::vector<int> > collectNeighbours(const std::vector<std::pair<int, int>>& IdPairs, std::vector<int> cellIDs);
	std::vector<int> getCellIDs(const std::vector<std::pair<int, int>>& IdPairs);
	std::vector<std::pair<int, int>> readVoronoiIndices(const std::string& filename);

	void readSnapshot(const std::string& snapshot);
	double squaredDistance(const std::vector<double>& point1, const std::vector<double>& point2);

protected:
	std::string fileMeshIndices;
	std::string snapshot;

	std::vector<std::pair<int, int>> IdPairs;

};

#endif /* MESH_H_ */
