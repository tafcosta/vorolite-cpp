#include "common_includes.h"
#include "Mesh.h"
#include "Rays.h"

void parseRayParamFile(const std::string& fileName, int& numRays, double& maxRadius,
                       std::vector<double>& sourceLocation, std::string& meshFile, std::string& snapFile);



int main() {

    int nRays = 0;
    double maxRadius = 0.0;
    std::vector<double> sourceLocation(3, 0.5);
    std::string meshFile, snapFile;

    parseRayParamFile("rayParam.txt", nRays, maxRadius, sourceLocation, meshFile, snapFile);

    if (nRays == 0 || maxRadius == 0.0 || meshFile.empty() || snapFile.empty()) {
        std::cerr << "Error: Missing or invalid parameters in rayParam.txt" << std::endl;
        return 1;
    }


	std::cout << "Starting VoroLite++ (Version 0.1)!" << std::endl;
	std::cout << "We are using " << nRays << " rays." << std::endl;

	Mesh *mesh = new Mesh("./output/tess_001_indices.dat", "./output/snap_001.hdf5");
	Rays *rays = new Rays(nRays, maxRadius, {0.5, 0.5, 0.5}, *mesh);

	rays->doRayTracing();
	rays->outputResults();
  
	delete mesh;
	delete rays;
	return 0;
}

void parseRayParamFile(const std::string& fileName, int& numRays, double& maxRadius,
                       std::vector<double>& sourceLocation, std::string& meshFile, std::string& snapFile) {
    std::ifstream inputFile(fileName);
    std::string line;

    if (!inputFile.is_open()) {
        std::cerr << "Error opening parameter file: " << fileName << std::endl;
        return;
    }

    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string key;
        std::string value;

        if (line.empty() || line[0] == '#') continue;

        std::getline(ss, key, '=');
        std::getline(ss, value);

        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);

        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);

        if (key == "numRays") {
            numRays = std::stoi(value);
        }
        else if (key == "maxRadius") {
            maxRadius = std::stod(value);
        }
        else if (key == "sourceLocation") {
            std::stringstream locStream(value);
            double x, y, z;
            char comma;
            locStream >> x >> comma >> y >> comma >> z;
            sourceLocation = {x, y, z};
        }
        else if (key == "meshFile") {
            meshFile = value;
        }
        else if (key == "snapFile") {
            snapFile = value;
        }
    }
}
