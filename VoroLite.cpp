#include "common_includes.h"
#include "Mesh.h"
#include "Rays.h"

void parseRayParamFile(const std::string& fileName, int& numRays, double& maxRadius,
                       std::vector<double>& sourceLocation, int& flowFilter, double& maxColumn, std::string& meshFile,
                       std::string& snapFile, std::string& ofileName);


// Parts of the main() function, as well as corresponding class functions like Rays.outputResults(), have been modified by L. Tortora
// Main changes: made the parameter file an input, to vary the inputs and outputs of VoroLite++
int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameter_file>" << std::endl;
        return 1;
    }

    std::string paramFile = argv[1];
    
    std::cout << "We are getting our parameters from \'" << paramFile << "\'" <<  std::endl;

    int flowFilter = 0;
    int nRays = 0;
    double maxRadius = 0.0;
    double maxColumn = 0.0;
    std::vector<double> sourceLocation(3, 0.5);
    std::string meshFile, snapFile, ofileName;

    parseRayParamFile(paramFile, nRays, maxRadius, sourceLocation, flowFilter, maxColumn, meshFile, snapFile, ofileName);

    if (nRays == 0 || maxRadius == 0.0 || maxColumn == 0.0 || meshFile.empty() || snapFile.empty()) {
        std::cerr << "Error: Missing or invalid parameters in rayParam.txt" << std::endl;
        return 1;
    }

	std::cout << "Starting VoroLite++ (Version 0.1)!" << std::endl;
	std::cout << "We are using " << nRays << " rays." << std::endl;
    std::cout << "The source is at position " << sourceLocation[0] << ", " << sourceLocation[1] << ", " << sourceLocation[2] << " (code units)" << std::endl;
    std::cout << "The maximum radius is " << maxRadius << " (code units)" << std::endl;

    Mesh *mesh = new Mesh(meshFile, snapFile);

    Rays *rays = new Rays(nRays, maxRadius, sourceLocation, flowFilter, maxColumn, *mesh);

	rays->doRayTracing();
	rays->outputResults(ofileName);
  
	delete mesh;
	delete rays;
	return 0;
}

void parseRayParamFile(const std::string& fileName, int& numRays, double& maxRadius,
                       std::vector<double>& sourceLocation, int& flowFilter, double& maxColumn, std::string& meshFile,
                       std::string& snapFile, std::string& ofileName) {
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
        else if (key == "flowFilter") {
        	flowFilter = std::stoi(value);
        }
        else if (key == "maxColumn") {
        	maxColumn = std::stod(value);
        }
        else if (key == "meshFile") {
            meshFile = value;
        }
        else if (key == "snapFile") {
            snapFile = value;
        }
        else if (key == "outputFile") {
            ofileName = value;
        }
    }
}
