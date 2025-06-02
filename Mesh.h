/*
 * Mesh.h
 *
 *  Created on: 19 Jan 2025
 *      Author: Tiago Costa (Newcastle University)
 */

#ifndef MESH_H_
#define MESH_H_

#include "common_includes.h"

class Mesh {
public:
	Mesh(std::string fileMeshIndices, std::string snapshot, double maxRadius, std::vector<double> sourcePosition);
	virtual ~Mesh();

	const double protonMass = 1.673e-24;

	int numCells;
	double boxSize;
	double maxRadius;
    std::vector<double> sourcePosition;

    double unitLength;
    double unitMass;
    double unitVelocity;

    std::vector<int> cellVisitsByRay;

    std::vector<double> cellFlux;
    std::vector<double> cellIncomingFlux;

    std::vector<double> cellLocalColumn;
    std::vector<double> cellHIIFraction;
    std::vector<std::vector<float>> cellCoordinates;
    std::vector<std::vector<float>> cellVelocities;

	std::vector<std::vector<int>> neighbourList;
	std::vector<int> findHostCellID(const std::vector<double>& target, int cellGuess);
	bool checkIfExitCellNeighboursCurrentCell(int iCell, int exitCell);

	double getDistanceToCell(const std::vector<double>& target, int cellIndex);
	double getDistanceBetweenCells(int iCell, int jCell);

	double getFluxOfRayInCell(int iRay, int iCell);
	void setFluxOfRayInCell(int iRay, int iCell, double newValue);

	double getHIIFraction(int iCell);
	double getMass(int iCell);
	double getDensity(int iCell);
	double getNumberDensity_in_cgs(int iCell);
	double getElectronNumberDensity_in_cgs(int iCell);
	double getMeanMolecularWeight(int iCell);

	int getIndex(int iCell);

	void setHIIFraction(int iCell, double newValue);
    void resetFluxes();
    void resizeFluxOfRayInCell(int iRay, int numVisitedCells);

private:
	std::vector<int> cellIndices;
	std::vector<int> cellIDs;
	
    std::vector<double> cellMass;
    std::vector<double> cellDensity;

    std::vector<std::vector<double>> fluxOfRayInCell;

	std::vector<std::vector<int> > collectNeighbours(const std::vector<std::pair<int, int>>& IdPairs, std::vector<int>& cellIDs);
	std::vector<int> getCellIDs(const std::vector<std::pair<int, int>>& IdPairs);
	std::vector<std::pair<int, int>> readVoronoiIndices(const std::string& filename);

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
