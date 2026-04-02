# VoroLite++

VoroLite++ is a C++ tool designed for ray tracing on a Voronoi grid. It performs ray tracing calculations and returns valuable information, including column densities and lists of densities along the rays. The tool works with Arepo snapshot files and requires a parameter file to configure the ray tracing process.
  
## Requirements

- C++17 or later.
- HDF5 library (for reading snapshot data).

## File Structure

- **rayParam.txt**: Parameter file that configures the ray tracing process.
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
   - `g++ -std=c++20 VoroLite.cpp Mesh.cpp Rays.cpp -lhdf5_cpp -lhdf5 -o vorolite`

4. The executable `vorolite` will be created.


## Usage

1. **Create a `rayParam.txt` file**: This file contains the parameters required for ray tracing. It should look like the following:

```txt
cosmo = false
nOutputs = 100
HIionisationCrossSection = 6.3e-18
HeIionisationCrossSection = 0.0
HeIIionisationCrossSection = 0.0
maxRadius = 0.12
sourceLocation = 0.5, 0.5, 0.5
lumTotal = 1.e54
timeMax  = 2e-6
dtime    = 5.e-10
meshFile = /Users/ntc132/eclipse-workspace/VoroLite++/output/tess_005_indices.dat
snapFile = /Users/ntc132/eclipse-workspace/VoroLite++/output/snap_005.hdf5
outputDirectory = /Users/ntc132/eclipse-workspace/VoroLite++/outputLight/
```

- `cosmo`: Set to false or 0 if simulation is not cosmological, or true or 1 for a cosmological simulation.
- `numOutputs`: Number of outputs.
- `HIionisationCrossSection`: Ionisation cross-section for neutral hydrogen.
- `HeIionisationCrossSection`: Ionisation cross-section for neutral helium.
- `HeIIionisationCrossSection`: Ionisation cross-section for singly ionised helium.
- `maxRadius`: Maximum radius for ray tracing.
- `sourceLocation`: The starting point of the rays (in x, y, z coordinates).
- `lumTotal`: Photon injection rate.
- `timeMax`: Final time in code units.
- `dtime`: Time step
- `meshFile`: Path to the Voronoi mesh file (typically `.dat`).
- `snapFile`: Path to the snapshot file (typically `.hdf5`).
- `outputDirectory`: Path to the output file.
- `initHIIFile` : (optional) Path to HII initialisation file,

2. **Run the ray tracing**:
- `./vorolite rayParam.txt`

3. **Output**: The program will generate an output file based on the one given as an input parameter in `rays_param.txt` (`outputFile`), which contains various results, including column densities, line of sight velocity along rays (weighted by column density), total distance travelled by rays, ray directions and number of cells traversed by each ray.

   
