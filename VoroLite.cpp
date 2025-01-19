#include "common_includes.h"
#include "Mesh.h"
#include "Rays.h"

int main() {

	std::cout << "Starting VoroLite++ (Version 0.1)!" << std::endl;

	Mesh *mesh = new Mesh("./output/tess_001_indices.dat");
	Rays *rays = new Rays(100, {0.5, 0.5, 0.5});

	std::cout << "All done!" << std::endl;
  
	delete mesh;
	return 0;
}
