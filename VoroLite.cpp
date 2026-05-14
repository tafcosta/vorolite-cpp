#include "common_includes.h"
#include "Mesh.h"
#include "Rays.h"
#include "Source.h"
#include "SourceVariable.h"

void parseRayParamFile(const std::string& fileName, bool& cosmo,
		int& nOutputs,
		double& HIionisationXsection,
		double& HeIionisationXsection,
		double& HeIIionisationXsection,
		double& dustAbsorptionOpacity,  double& maxRadius,
        int64_t& Nside, std::string& meshFile,
        std::string& snapFile, std::string& oDirectory);

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameter_file>" << std::endl;
        return 1;
    }

    std::string paramFile = argv[1];
    std::cout << "We are getting our parameters from \'" << paramFile << "\'" <<  std::endl;

    int nOutputs = 100;

    double HIionisationCrossSection   = 0.0;
    double HeIionisationCrossSection  = 0.0;
    double HeIIionisationCrossSection = 0.0;
    double dustAbsorptionOpacity      = 0.0;

    double maxRadius = 0.0;

    int64_t Nside = 4;
    bool cosmo = false;
    std::string meshFile, snapFile, oDirectory;

    parseRayParamFile(paramFile, cosmo, nOutputs, HIionisationCrossSection, HeIionisationCrossSection, HeIIionisationCrossSection, dustAbsorptionOpacity,
    		maxRadius, Nside, meshFile, snapFile, oDirectory);

    if (maxRadius == 0.0 || meshFile.empty() || snapFile.empty()) {
        std::cerr << "Error: Missing or invalid parameters in rayParam.txt" << std::endl;
        return 1;
    }

	std::cout << "Starting HeatThatDust!" << std::endl;
    Mesh *mesh = new Mesh(meshFile, snapFile, maxRadius, cosmo);
    Rays *rays = new Rays(HIionisationCrossSection, HeIionisationCrossSection, HeIIionisationCrossSection, maxRadius, Nside, *mesh);

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

    int snapshotIndex = 0;
    std::cout << "Setting up rays..." << std::endl;
    rays->calculateRays();
    std::cout << "Setting up rays OK." << std::endl;

    std::ostringstream ofName;
    ofName << oDirectory << "rays_output_" << snapshotIndex << ".txt";
    std::string ofileName = ofName.str();
    rays->outputResults(ofileName);

    std::cout << "Starting radiative transfer" << std::endl;
    mesh->resetPhotons();
    rays->doRadiativeTransfer();

    if (outFile.is_open()) {

    	outFile << std::setprecision(15) << std::scientific;
    	for (int iCell = 0; iCell < mesh->numCells; ++iCell) {
    		outFile << mesh->getIndex(iCell) << " ";
    		for (float coord : mesh->cellCoordinates[iCell]) {
    			outFile << coord << " ";
    		}
    		outFile << mesh->getHIIFraction(iCell) << " "
    				<< mesh->getHeIIFraction(iCell) << " "
					<< mesh->getHeIIIFraction(iCell) << " "
					<< mesh->cellIncomingPhotonRate[iCell]
					<< std::endl;
                }
    } else {
    	std::cerr << "Unable to open file " << filename.str() << " for writing." << std::endl;
    }

	delete mesh;
	delete rays;

	return 0;
}

void parseRayParamFile(const std::string& fileName, bool& cosmo, int& nOutputs,
		double& HIionisationCrossSection, double& HeIionisationCrossSection, double& HeIIionisationCrossSection,
		double& dustAbsorptionOpacity, double& maxRadius,
        int64_t& Nside, std::string& meshFile,
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

        if (key == "cosmo") {
            if (value == "true" || value == "1") {
                cosmo = true;
            }
            else if (value == "false" || value == "0") {
                cosmo = false;
            }
            else {
                throw std::runtime_error("Invalid value for cosmo: " + value);
            }
        }
        else if (key == "nOutputs") {
        	nOutputs = std::stoi(value);
        }
        else if (key == "HIionisationCrossSection") {
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
