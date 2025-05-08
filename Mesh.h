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
	Mesh(std::string fileMeshIndices, std::string snapshot, double maxRadius, std::vector<double> sourceLocation);
	virtual ~Mesh();

	int numCells;
	double boxSize;
	double maxRadius;


	std::vector<int> cellIDs;
    std::vector<std::vector<float>> cellCoordinates;
    std::vector<std::vector<float>> cellVelocities;
    std::vector<double> cellDensity;
    std::vector<double> cellFlux;

	std::vector<std::vector<int>> neighbourList;
	std::vector<int> findHostCellID(const std::vector<double>& target, int cellGuess);
	bool checkIfExitCellNeighboursCurrentCell(int iCell, int exitCell);

	double getDistanceToCell(const std::vector<double>& target, int cellIndex);
	double getDistanceBetweenCells(int iCell, int jCell);


private:
	std::vector<std::vector<int> > collectNeighbours(const std::vector<std::pair<int, int>>& IdPairs, std::vector<int>& cellIDs);
	std::vector<int> getCellIDs(const std::vector<std::pair<int, int>>& IdPairs);
	std::vector<std::pair<int, int>> readVoronoiIndices(const std::string& filename);

    std::vector<double> sourcePosition;


	void getNumCellsInRegion();
	void saveVoronoiIndices(const std::string& filename, const std::vector<std::pair<int, int>>& IdPairs);
	void readSnapshot(const std::string& snapshot);
	double squaredDistance(const std::vector<float>& point1, const std::vector<double>& point2);

protected:
	std::string fileMeshIndices;
	std::string snapshot;

	std::vector<std::pair<int, int>> IdPairs;

};

#endif /* MESH_H_ */
