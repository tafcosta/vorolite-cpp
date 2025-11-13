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
	double HubbleParam;
	double scaleFactor;

	double xHydrogen = 0.76;
	double xHelium   = 0.24;

    std::vector<int> cellVisitsByRay;

    std::vector<double> cellFlux;
    std::vector<double> cellIncomingFlux;

    std::vector<double> cellLocalColumn;
    std::vector<double> cellHIIFraction;
	std::vector<double> cellHIFraction;

    std::vector<double> cellHeIIFraction;
    std::vector<double> cellHeIIIFraction;

	std::vector<double> cellElectronFraction;
	std::vector<double> cellXH;
	std::vector<double> cellMetallicity;
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
	double getHeIIFraction(int iCell);
	double getHeIIIFraction(int iCell);

	double getMass(int iCell);
	double getDensity(int iCell);
	double getHNumberDensity_in_cgs(int iCell);
	double getHeNumberDensity_in_cgs(int iCell);

	double getElectronNumberDensity_in_cgs(int iCell);
	double getMeanMolecularWeight(int iCell);
	double getMetallicityInSolar(int iCell);
	double getSelfShieldingCorrection(int iCell);
    double getFlux(int iCell);
    double getIncomingFlux(int iCell);

	int getIndex(int iCell);

	void setHIIFraction(int iCell, double newValue);
	void setHeIIFraction(int iCell, double newValue);
	void setHeIIIFraction(int iCell, double newValue);

    void resetFluxes();
    void resizeFluxOfRayInCell(int iRay, int numVisitedCells);
	void doSelfShieldingCorrection();

private:
	std::vector<int> cellIndices;
	std::vector<int> cellIDs;
    std::vector<double> cellMass;
    std::vector<double> cellDensity;

    std::vector<std::vector<double>> fluxOfRayInCell;

    void appendDensity(H5::H5File& file);
    void appendMass(H5::H5File& file);
    void appendIDs(H5::H5File& file);
    void appendCoordinates(H5::H5File& file);
    void appendVelocities(H5::H5File& file);
	void appendHIFraction(H5::H5File& file);
	void appendElectronFraction(H5::H5File& file);
	void appendXH(H5::H5File& file);
    void appendMetallicity(H5::H5File& file);
    void readHeader(H5::H5File& file);
    std::vector<std::string> getSnapshotFiles(const std::string& snapshotBase);

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
