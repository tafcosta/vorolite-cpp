# VoroLite++

VoroLite++ is a C++ tool designed for ray tracing on a Voronoi grid. It performs ray tracing calculations and returns valuable information, including column densities and lists of densities along the rays. The tool works with snapshot files and requires a parameter file to configure the ray tracing process.
  
## Requirements

- C++17 or later.
- HDF5 library (for reading snapshot data).

## File Structure

- **rays_param.txt**: Parameter file that configures the ray tracing process.
- **meshFile**: Mesh file containing the Voronoi tessellation indices.
- **snapFile**: Snapshot file containing the data for ray tracing.
- **rays_output.txt**: Output file containing list of rays, ray directions and column densities.


## Installation

To compile and run VoroLite++, follow these steps:

1. Clone the repository:
   - `git clone https://github.com/your_username/vorolite_cpp.git`
   - `cd vorolite_cpp`

2. Ensure you have a C++20 compatible compiler and the required libraries (HDF5) installed.

3. Compile the code:
   - `g++ -std=c++20 -lhdf5_cpp -lhdf5  -o vorolite vorolite++.cpp`

4. The executable `vorolite` will be created.


## Usage

1. **Create a `rays_param.txt` file**: This file contains the parameters required for ray tracing. It should look like the following:

```txt
numRays = 10000
maxRadius = 0.5
sourceLocation = 0.5, 0.5, 0.5
meshFile = ./output/tess_001_indices.dat
snapFile = ./output/snap_001.hdf5
```

- `numRays`: Number of rays to trace.
- `maxRadius`: Maximum radius for ray tracing.
- `sourceLocation`: The starting point of the rays (in x, y, z coordinates).
- `meshFile`: Path to the Voronoi mesh file (typically `.dat`).
- `snapFile`: Path to the snapshot file (typically `.hdf5`).

2. **Run the ray tracing**:
- `./vorolite rays_param.txt`

3. **Output**: The program will generate an output file named `ray_output.txt`, which contains various results, including column densities and a list of densities along each ray.

   
