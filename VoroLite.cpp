#include "common_includes.h"
#include "Mesh.h"

int main() {

	std::cout << "Starting VoroLite++ (Version 0.1)!" << std::endl;

	Mesh *mesh = new Mesh("./output/tess_001_indices.dat");

	std::cout << "All done!" << std::endl;
  
	delete mesh;
	return 0;
}
