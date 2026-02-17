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
    cellPhotonRate.resize(numCells, 0.0);
    cellIncomingPhotonRate.resize(numCells, 0.0);

    cellAbsorbedPhotonRateHI.resize(numCells, 0.0);
    cellAbsorbedPhotonRateHeI.resize(numCells, 0.0);
    cellAbsorbedPhotonRateHeII.resize(numCells, 0.0);

    cellNetIonisationRate.resize(numCells, 0.0);

	cellHIIFraction.resize(numCells, 0.0);
	cellHeIIFraction.resize(numCells, 0.0);
	cellHeIIIFraction.resize(numCells, 0.0);

	cellRemainingHI.resize(numCells, 0.0);
	cellRemainingHeI.resize(numCells, 0.0);
	cellRemainingHeII.resize(numCells, 0.0);

	fluxOfRayInCell.resize(numCells);       //The first dimension should be number of rays

    //doSelfShieldingCorrection();

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

double Mesh::getSpecificInternalEnergy(int iCell){
	return cellSpecificInternalEnergy[iCell];
}

double Mesh::getHNumberDensity_in_cgs(int iCell){
	return xHydrogen * cellDensity[iCell] / protonMass * (unitMass / (scaleFactor * unitLength * scaleFactor * unitLength * scaleFactor * unitLength)) * HubbleParam * HubbleParam;
}

double Mesh::getHeNumberDensity_in_cgs(int iCell) {
    double heliumMass = 4.0 * protonMass;
    return yHelium * cellDensity[iCell] / heliumMass * (unitMass / (scaleFactor * unitLength * scaleFactor * unitLength * scaleFactor * unitLength)) * HubbleParam * HubbleParam;
}

double Mesh::getElectronNumberDensity_in_cgs(int iCell){
	return cellHIIFraction[iCell] * getHNumberDensity_in_cgs(iCell);
}

double Mesh::getMeanMolecularWeight(int iCell){
    double xHII   = getHIIFraction(iCell);
    double yHeII  = getHeIIFraction(iCell);
    double zHeIII = getHeIIIFraction(iCell);

    double muInv = xHydrogen * (1.0 + xHII) + (yHelium / 4.0) * (1.0 + yHeII + 2.0 * zHeIII);
    return 1.0 / muInv;
}

double Mesh::getTemperature_in_K(int iCell){
	return 10000.;/*getSpecificInternalEnergy(iCell) * unitVelocity * unitVelocity *
			(adiabaticIndex - 1.0) * getMeanMolecularWeight(iCell) * protonMass / boltzmannConstant;*/
}

double Mesh::getMetallicityInSolar(int iCell){
	return cellMetallicity[iCell]/0.0127;
}

double Mesh::getRemainingHI(int iCell){
	return cellRemainingHI[iCell];
}

double Mesh::getRemainingHeI(int iCell){
	return cellRemainingHeI[iCell];
}

double Mesh::getRemainingHeII(int iCell){
	return cellRemainingHeII[iCell];
}

double Mesh::getSelfShieldingCorrection(int iCell) {
    const double rho_s = 1.52e-2;
    const double rho_u = 4.53e-3;
    const double p     = 2.68;

    double nH   = cellXH[iCell] * getHNumberDensity_in_cgs(iCell);
    double f_hi = 1 - cellHIIFraction[iCell];

    double new_f_hi = f_hi;

    if (nH >= rho_u && nH <= rho_s) {
        double numerator   = f_hi * std::pow(rho_s - nH, p)
                           + std::pow(nH - rho_u, p);
        double denominator = std::pow(rho_s - rho_u, p);
        new_f_hi = numerator / denominator;
    } 
    else if (nH > rho_s) {
        new_f_hi = 1.0;
    }

    return new_f_hi;
}

void Mesh::doSelfShieldingCorrection() {
    for (int iCell = 0; iCell < numCells; ++iCell) {
        double newcellHIFraction = getSelfShieldingCorrection(iCell);
        cellHIIFraction[iCell] = 1.0 - newcellHIFraction;
    }
}

double Mesh::getFluxOfRayInCell(int iRay, int iCell){
	return fluxOfRayInCell[iRay][iCell];
}

double Mesh::getIncomingPhotonRate(int iCell){
	return cellIncomingPhotonRate[iCell];
}

double Mesh::getAbsorbedPhotonRateHI(int iCell){
	return cellAbsorbedPhotonRateHI[iCell];
}

double Mesh::getHIIFraction(int iCell){
	return cellHIIFraction[iCell];
}

double Mesh::getHeIIFraction(int iCell){
	return cellHeIIFraction[iCell];
}

double Mesh::getHeIIIFraction(int iCell){
	return cellHeIIIFraction[iCell];
}

int Mesh::getIndex(int iCell){
	return cellIndices[iCell];
}

void Mesh::setFluxOfRayInCell(int iRay, int iCell, double newValue){
	fluxOfRayInCell[iRay][iCell] = newValue;
}

void Mesh::setHIIFraction(int iCell, double newValue){
	cellHIIFraction[iCell] = newValue;

	if(newValue > 1.0)
		cellHIIFraction[iCell] = 1.0;

	if(newValue < 0.0)
		cellHIIFraction[iCell] = 0.0;
}

void Mesh::setHeIIFraction(int iCell, double newValue){
	cellHeIIFraction[iCell] = newValue;

	if(newValue > 1.0)
		cellHeIIFraction[iCell] = 1.0;

	if(newValue < 0.0)
		cellHeIIFraction[iCell] = 0.0;
}

void Mesh::setHeIIIFraction(int iCell, double newValue){
	cellHeIIIFraction[iCell] = newValue;

	if(newValue > 1.0)
		cellHeIIIFraction[iCell] = 1.0;

	if(newValue < 0.0)
		cellHeIIIFraction[iCell] = 0.0;
}

void Mesh::setRemainingHI(int iCell, double newValue){
	cellRemainingHI[iCell] = std::max(0.0, newValue);
}

void Mesh::setRemainingHeI(int iCell, double newValue){
	cellRemainingHeI[iCell] = std::max(0.0, newValue);
}

void Mesh::setRemainingHeII(int iCell, double newValue){
	cellRemainingHeII[iCell] = std::max(0.0, newValue);
}

void Mesh::getNumCellsInRegion(){
    std::vector<std::vector<float>> filteredCoordinates;
    std::vector<std::vector<float>> filteredVelocities;
    std::vector<double> filteredDensity;
    std::vector<double> filteredSpecificInternalEnergy;
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
            filteredSpecificInternalEnergy.push_back(cellSpecificInternalEnergy[iCell]);
            filteredMasses.push_back(cellMass[iCell]);
            filteredIDs.push_back(cellIDs[iCell]);
            filteredCellIndices.push_back(cellIndices[iCell]);
		}
	}

    cellCoordinates = std::move(filteredCoordinates);
    cellVelocities = std::move(filteredVelocities);
    cellDensity = std::move(filteredDensity);
    cellSpecificInternalEnergy = std::move(filteredSpecificInternalEnergy);
    cellIDs = std::move(filteredIDs);
    cellIndices = std::move(filteredCellIndices);
    cellMass = std::move(filteredMasses);

    numCells = cellDensity.size();

    std::cout << "Reduced to " << numCells << " cells within maxRadius = " << maxRadius << std::endl;
}

void Mesh::resetPhotons(){
    for(int iCell = 0; iCell < numCells; iCell++){
        cellPhotonRate[iCell] = 0.;

        cellAbsorbedPhotonRateHI[iCell]   = 0.;
        cellAbsorbedPhotonRateHeI[iCell]  = 0.;
        cellAbsorbedPhotonRateHeII[iCell] = 0.;

        cellIncomingPhotonRate[iCell] = 0.;
        cellNetIonisationRate[iCell]  = 0.;

        const double xMax = 1.0 - 1e-20;

        double xHII = getHIIFraction(iCell);
        if (xHII < 0.0) xHII = 0.0;
        if (xHII > xMax) xHII = xMax;

        double yHeII = getHeIIFraction(iCell);
        double zHeIII = getHeIIIFraction(iCell);
        if (yHeII < 0.0) yHeII = 0.0;
        if (zHeIII < 0.0) zHeIII = 0.0;
        if (yHeII + zHeIII > xMax) {
            const double s = yHeII + zHeIII;
            yHeII = yHeII * xMax / s;
            zHeIII = zHeIII * xMax / s;
        }

        const double neutralH = std::max(1.0 - xHII, 1e-20);
        const double neutralHe = std::max(1.0 - yHeII - zHeIII, 1e-20);

        double volume = getMass(iCell) / getDensity(iCell) *
             scaleFactor * unitLength  *
             scaleFactor * unitLength  *
             scaleFactor * unitLength  *
             HubbleParam * HubbleParam * HubbleParam;

        double nH  = getHNumberDensity_in_cgs(iCell);
        double nHe = getHeNumberDensity_in_cgs(iCell);

        cellRemainingHI[iCell]  = neutralH * (nH * volume);
        cellRemainingHeI[iCell] = neutralHe * (nHe * volume);
        cellRemainingHeII[iCell]= yHeII * (nHe * volume);
    }
}



void Mesh::readSnapshot(const std::string& snapshotBase) {
    try {
        // std::cout << "Getting into getSnapshotFiles" << std::endl;
    	// std::vector<std::string> files = getSnapshotFiles(snapshotBase);

    	// std::cout << "[DEBUG] getSnapshotFiles() returned " << files.size() << " files:" << std::endl;
    	// for (const auto& f : files) {
    	//     std::cout << " - " << f << std::endl;
    	// }

        // // std::vector<std::string> files = snapshotBase;

        // bool headerRead = false;

        // for (const std::string& fileName : files) {
        //     H5::H5File file(fileName, H5F_ACC_RDONLY);

        //     if (!headerRead) {
        //         readHeader(file);
        //         headerRead = true;
        //     }

        //     appendDensity(file);
        //     appendMass(file);
        //     appendIDs(file);
        //     appendCoordinates(file);
        //     appendVelocities(file);
        // }

        bool headerRead = false;
        bool cosmo = false;

        H5::H5File file(snapshotBase, H5F_ACC_RDONLY);
        if (!headerRead) {
            readHeader(file);
            headerRead = true;
        }

        if(!cosmo){
        	scaleFactor = 1.0;
        	HubbleParam = 1.0;
        }

        appendDensity(file);
        appendSpecificInternalEnergy(file);
        appendMass(file);
        appendIDs(file);
        appendCoordinates(file);
        appendVelocities(file);
        //appendMetallicity(file);
        //appendElectronFraction(file);
        //appendXH(file);

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
    headerGroup.openAttribute("HubbleParam").read(H5::PredType::NATIVE_DOUBLE, &HubbleParam);
    headerGroup.openAttribute("Time").read(H5::PredType::NATIVE_DOUBLE, &scaleFactor);
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

void Mesh::appendSpecificInternalEnergy(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/InternalEnergy");
    H5::DataSpace space = dataset.getSpace();

    hsize_t numElements;
    space.getSimpleExtentDims(&numElements);

    std::vector<double> buffer(numElements);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    cellSpecificInternalEnergy.insert(cellSpecificInternalEnergy.end(), buffer.begin(), buffer.end());
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

/*
void Mesh::appendMetallicity(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/GFM_Metallicity");
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t numElements;
    dataspace.getSimpleExtentDims(&numElements);

    size_t offset = cellMetallicity.size();
    cellMetallicity.resize(offset + numElements);

    dataset.read(cellMetallicity.data() + offset, H5::PredType::NATIVE_DOUBLE);
}
*/

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

/*
void Mesh::appendElectronFraction(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/ElectronAbundance");
    H5::DataSpace space = dataset.getSpace();

    hsize_t numElements;
    space.getSimpleExtentDims(&numElements);

    std::vector<double> buffer(numElements);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    cellElectronFraction.insert(cellElectronFraction.end(), buffer.begin(), buffer.end());
}
*/

/*
void Mesh::appendXH(H5::H5File& file) {
    H5::DataSet dataset = file.openDataSet("/PartType0/GFM_Metals");
    H5::DataSpace space = dataset.getSpace();

    hsize_t dims[2];
    space.getSimpleExtentDims(dims);

    std::vector<float> buffer(dims[0] * dims[1]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_FLOAT);

    for (hsize_t i = 0; i < dims[0]; ++i) {
        double XH = buffer[i * dims[1] + 0]; // Assuming Hydrogen is the first element
        cellXH.push_back(XH);
    }
}
*/

std::vector<std::string> Mesh::getSnapshotFiles(const std::string& snapshotPath) {
    std::vector<std::string> files;

    std::filesystem::path inputPath(snapshotPath);
    std::filesystem::path dir = inputPath.parent_path();
    if (dir.empty()) dir = ".";

    // std::cout << "Trying baseName" << std::endl;
    // std::cout << "baseName" << inputPath << std::endl;
    // std::string baseName = inputPath.stem().string();     // "snap_006"
    // std::cout << "Trying extension" << std::endl;
    // std::string extension = inputPath.extension().string(); // ".hdf5"

    // std::cout << "Trying step 1" << std::endl;
    // // Step 1: Look for split files like snap_006.0.hdf5, snap_006.1.hdf5, ...
    // for (const auto& entry : std::filesystem::directory_iterator(dir)) {
    //     if (!entry.is_regular_file()) continue;

    //     std::string fname = entry.path().filename().string();

    //     // Match baseName.N.hdf5, e.g., snap_006.0.hdf5
    //     if (fname.rfind(baseName + ".", 0) == 0 &&
    //         fname.size() > baseName.size() + 6 &&
    //         fname.compare(fname.size() - 6, 6, ".hdf5") == 0) {
    //         files.push_back(entry.path().string());
    //     }
    // }

    std::cout << "Going directly to step 2" << std::endl;
    std::cout << "Trying with " << snapshotPath << std::endl;
    std::cout << "Empty? " << files.empty() << std::endl;
    std::cout << "Exists? " << std::filesystem::exists(snapshotPath) << std::endl;
    // Step 2: Fallback to monolithic file if no split files found
    if (files.empty() && std::filesystem::exists(snapshotPath)) {
        std::cout << "Made it here" << std::endl;
        files.push_back(snapshotPath);
    }
    std::cout << "Got through step 2" << std::endl;

    std::sort(files.begin(), files.end());
    return files;
}

Mesh::~Mesh() {

}
