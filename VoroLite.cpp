#include "common_includes.h"
#include "Mesh.h"
#include "Photochemistry.h"
#include "Rays.h"
#include "Source.h"

void parseRayParamFile(const std::string& fileName, double& HIionisationXsection,
		double& HeIionisationXsection,
		double& HeIIionisationXsection,
		double& dustAbsorptionOpacity,  double& maxRadius,
        std::vector<double>& sourcePosition, double& lumTotal, double& timeMax, int64_t& Nside, std::string& meshFile,
        std::string& snapFile, std::string& oDirectory);

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameter_file>" << std::endl;
        return 1;
    }

    std::string paramFile = argv[1];
    std::cout << "We are getting our parameters from \'" << paramFile << "\'" <<  std::endl;

    double HIionisationCrossSection     = 0.0;
    double HeIionisationCrossSection    = 0.0;
    double HeIIionisationCrossSection   = 0.0;

    double dustAbsorptionOpacity = 0.0;

    double maxRadius = 0.0;
    double lumTotal  = 0.0;
    double timeMax   = 0.0;
    int64_t Nside = 4;
    std::vector<double> sourcePosition(3, 0.5);
    std::string meshFile, snapFile, oDirectory;
    std::filesystem::path lightcurvefile = "data/Lion_basic_ref.txt";

    parseRayParamFile(paramFile, HIionisationCrossSection, HeIionisationCrossSection, HeIIionisationCrossSection, dustAbsorptionOpacity, maxRadius, sourcePosition, lumTotal, timeMax, Nside, meshFile, snapFile, oDirectory);

    if (maxRadius == 0.0 || meshFile.empty() || snapFile.empty()) {
        std::cerr << "Error: Missing or invalid parameters in rayParam.txt" << std::endl;
        return 1;
    }

	std::cout << "Starting VoroLite++ RT (Version 0.1)!" << std::endl;

    std::cout << "Mesh initialisation starting..." << std::endl;
    Mesh *mesh = new Mesh(meshFile, snapFile, maxRadius, sourcePosition);
    std::cout << "Mesh initialisation OK." << std::endl;

    std::cout << "Source initialisation starting..." << std::endl;
    Source *source = new Source(sourcePosition, lumTotal);
    std::cout << "Source initialisation OK." << std::endl;

    std::cout << "Rays initialisation starting..." << std::endl;
    Rays *rays = new Rays(HIionisationCrossSection, maxRadius, sourcePosition, lumTotal, Nside, *mesh, *source);
    std::cout << "Rays initialisation OK." << std::endl;

    std::cout << "Photochemistry initialisation starting..." << std::endl;
    Photochemistry *photochemistry = new Photochemistry(*mesh, *rays, HIionisationCrossSection, HeIionisationCrossSection, HeIIionisationCrossSection);
    std::cout << "Photochemistry initialisation OK" << std::endl;

    std::ostringstream filename;
    filename << oDirectory << "HIIfraction_init.txt";

    std::ofstream outFile(filename.str());
    if (outFile.is_open()) {
        for (int iCell = 0; iCell < mesh->numCells; ++iCell) {
            outFile << mesh->getIndex(iCell) << " ";
            for (float coord : mesh->cellCoordinates[iCell]) {
                outFile << coord << " ";
            }
            outFile << mesh->getHIIFraction(iCell) << " " << mesh->cellIncomingPhotonRate[iCell] << std::endl;
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file " << filename.str() << " for writing." << std::endl;
    }

    double time = 0;
    double dtime  = 5e-10;

    double printInterval = timeMax/100;
    double TimeNextOutput = printInterval;

    int snapshotIndex = 0;

    std::cout << "Setting up rays..." << std::endl;
    rays->calculateRays();
    std::cout << "Setting up rays OK." << std::endl;


    std::ostringstream ofName;
    ofName << oDirectory << "rays_output_" << snapshotIndex << ".txt";
    std::string ofileName = ofName.str();
    rays->outputResults(ofileName);

    std::cout << "Starting radiative transfer" << std::endl;


    while (time < timeMax) {

        mesh->resetPhotons();
    	rays->doRadiativeTransfer(time, dtime);
    	photochemistry->evolveIonisation(dtime);


        if (time >= TimeNextOutput) {
            std::cout << "time = " << time << std::endl;

            std::ostringstream filename;
            filename << oDirectory << "HIIfraction_" << snapshotIndex << ".txt";

            std::ofstream outFile(filename.str());
            if (outFile.is_open()) {
                for (int iCell = 0; iCell < mesh->numCells; ++iCell) {
                	outFile << mesh->getIndex(iCell) << " ";
                    for (float coord : mesh->cellCoordinates[iCell]) {
                        outFile << coord << " ";
                    }
                    outFile << mesh->getHIIFraction(iCell) << " " << mesh->getHeIIFraction(iCell) << " " << mesh->getHeIIIFraction(iCell) << " " << mesh->cellIncomingPhotonRate[iCell] << std::endl;
                }
                outFile.close();
            } else {
                std::cerr << "Unable to open file " << filename.str() << " for writing." << std::endl;
            }

            ++snapshotIndex;
            TimeNextOutput += printInterval;
        }

        time += dtime;
    }
  
	delete mesh;
	delete rays;
	delete source;
	delete photochemistry;

	return 0;
}

void parseRayParamFile(const std::string& fileName, double& HIionisationCrossSection, double& HeIionisationCrossSection, double& HeIIionisationCrossSection,
		double& dustAbsorptionOpacity, double& maxRadius,
        std::vector<double>& sourceLocation, double& lumTotal, double& timeMax, int64_t& Nside, std::string& meshFile,
        std::string& snapFile, std::string& oDirectory) {

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

        if (key == "HIionisationCrossSection") {
        	HIionisationCrossSection = std::stod(value);
        }
        else if (key == "HeIionisationCrossSection") {
        	HeIionisationCrossSection = std::stod(value);
        }
        else if (key == "HeIIionisationCrossSection") {
        	HeIIionisationCrossSection = std::stod(value);
        }
        else if (key == "dustAbsorptionOpacity") {
        	dustAbsorptionOpacity = std::stod(value);
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
        else if (key == "lumTotal") {
            lumTotal = std::stod(value);
        }
        else if (key == "timeMax") {
            timeMax = std::stod(value);
        }
        else if (key == "meshFile") {
            meshFile = value;
        }
        else if (key == "snapFile") {
            snapFile = value;
        }
        else if (key == "outputDirectory") {
            oDirectory = value;
        }
        else if (key == "Nside") {
            Nside = std::stoll(value);
        }
    }
}
