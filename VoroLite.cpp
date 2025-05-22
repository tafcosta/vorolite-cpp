#include "common_includes.h"
#include "Mesh.h"
#include "Photochemistry.h"
#include "Rays.h"

void parseRayParamFile(const std::string& fileName, double &Xsection, double& maxRadius,
                       std::vector<double>& sourceLocation, std::string& meshFile,
                       std::string& snapFile, std::string& ofileName);


int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameter_file>" << std::endl;
        return 1;
    }

    std::string paramFile = argv[1];
    
    std::cout << "We are getting our parameters from \'" << paramFile << "\'" <<  std::endl;

    double crossSection = 0.0;
    double maxRadius = 0.0;
    std::vector<double> sourcePosition(3, 0.5);
    std::string meshFile, snapFile, ofileName;

    parseRayParamFile(paramFile, crossSection, maxRadius, sourcePosition, meshFile, snapFile, ofileName);

    if (maxRadius == 0.0 || meshFile.empty() || snapFile.empty()) {
        std::cerr << "Error: Missing or invalid parameters in rayParam.txt" << std::endl;
        return 1;
    }

	std::cout << "Starting VoroLite++ RT (Version 0.1)!" << std::endl;

    Mesh *mesh = new Mesh(meshFile, snapFile, maxRadius, sourcePosition);
    Rays *rays = new Rays(crossSection, maxRadius, sourcePosition, *mesh);
    Photochemistry *photochemistry = new Photochemistry(*mesh, crossSection);

    std::cout << "The maximum radius is " << maxRadius << " (code units)" << std::endl;


    double time = 0;
    double timeMax = 0.001;
    double dtime = 0.0001;

    double printInterval = timeMax/5;
    double nextPrintTime = printInterval;

    while(time < timeMax){
    	rays->doRayTracing();
    	photochemistry->evolveIonisation(dtime);

        if (time >= nextPrintTime) {
            std::cout << "time = " << time << std::endl;
            nextPrintTime += printInterval;
        }

    	time += dtime;
    }


    std::ofstream outFile("HIIfraction.txt");
    if (outFile.is_open()) {
        for (int iCell = 0; iCell < mesh->numCells; ++iCell) {
            for (float coord : mesh->cellCoordinates[iCell]) {
                outFile << coord << " ";
            }
            outFile << mesh->cellHIIFraction[iCell] << " " << mesh->cellFlux[iCell] << std::endl;
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
  
	delete mesh;
	delete rays;
	delete photochemistry;

	return 0;
}

void parseRayParamFile(const std::string& fileName, double& crossSection, double& maxRadius,
                       std::vector<double>& sourceLocation, std::string& meshFile,
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

        if (key == "crossSection") {
        	crossSection = std::stod(value);
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
        else if (key == "outputFile") {
            ofileName = value;
        }
    }
}
