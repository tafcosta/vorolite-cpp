#include "common_includes.h"
#include "Mesh.h"
#include "Rays.h"

int main() {

	std::cout << "Starting VoroLite++ (Version 0.1)!" << std::endl;

	Mesh *mesh = new Mesh("./output/tess_001_indices.dat", "./output/snap_001.hdf5");
	Rays *rays = new Rays(1, {0.5, 0.5, 0.5}, *mesh);

	rays->doRayTracing();

	std::cout << "All done!" << std::endl;
  
	delete mesh;
	delete rays;
	return 0;
}
