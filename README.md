# BloodFlow Simulation

**BloodFlow** is a C++ simulation project that models blood flow dynamics at the single-cell level. It leverages the power of **Palabos** for lattice Boltzmann simulations, **LAMMPS** for molecular dynamics, and **Ascent** for advanced visualization. This project integrates these tools to provide a comprehensive framework for simulating and analyzing blood flow in biological systems.

## Features

- **Palabos Integration**: Utilizes Palabos for high-performance lattice Boltzmann simulations of incompressible fluid flows.
- **LAMMPS Coupling**: Integrates LAMMPS for molecular dynamics, enabling detailed modeling of cellular interactions and forces.
- **Ascent Visualization**: Employs Ascent for real-time visualization and analysis of simulation data, facilitating intuitive understanding of blood flow dynamics.
- **MPI Support**: Supports parallel computing using MPI to leverage multi-core and distributed computing environments.
- **Customizable Parameters**: Offers flexible configuration options to tailor simulations to specific research needs.

## Prerequisites

Before building and running the BloodFlow simulation, ensure that the following dependencies are installed on your system:

- **CMake**: Version 3.18 or higher
- **C++ Compiler**: Supporting C++11 standard (e.g., GCC, Clang, MSVC)
- **MPI**: For parallel execution (e.g., OpenMPI, MPICH)
- **Ascent**: For visualization (optional but recommended)
- **Git**: For cloning the repository and fetching dependencies

## Installation

### Clone the Repository

Begin by cloning the BloodFlow repository from GitLab (replace `<repository_url>` with the actual URL if different):

```bash
git clone <repository_url> BloodFlow
cd BloodFlow
```

### Configure the Build
BloodFlow uses CMake for build configuration. You can choose to let CMake automatically fetch and manage dependencies like Palabos and LAMMPS or use manually installed versions.

#### Options

- `FETCH_PALABOS` (ON by default): Download and manage Palabos via FetchContent.
- `FETCH_LAMMPS` (ON by default): Download and manage LAMMPS via ExternalProject.
- `BUILD_EXAMPLE_SINGLECELL` (ON by default): Build the SingleCell example.
- `ENABLE_ASCENT` (ON by default): Enable Ascent for visualization.
- `ENABLE_MPI` (ON by default): Enable MPI support.
- `ENABLE_POSIX` (ON by default on non-Windows systems): Enable POSIX support.

#### Example Configuration

Create a build directory and configure the project with default options:
```bash
mkdir build
cd build
cmake -DAscent_DIR=/path/to/lib/cmake/ascent ..
```

**Custom Configuration**: To customize build options, use `-D` flags. For example, to disable fetching LAMMPS and Ascent:
```bash
cmake -DFETCH_LAMMPS=OFF -DENABLE_ASCENT=OFF ..
```

**Specifying Installed Dependencies**: If you have Palabos or LAMMPS installed manually, disable fetching and set their paths:
```bash
cmake -DFETCH_PALABOS=OFF -DPALABOS_ROOT=/path/to/palabos \
      -DFETCH_LAMMPS=OFF -DLAMMPS_DIR=/path/to/lammps ..
```

#### Build the Project
Once configured, build the project using `make` (or your preferred build system):
```bash
cmake --build . --parallel
```
This will compile the BloodFlow simulation along with its dependencies if fetching is enabled.

### Running the Simulation
After a successful build, navigate to the `examples/singleCell` directory and execute the simulation:
```bash
cd examples/singleCell
mpirun -np <num_processes> ./singleCell <lammps_input_file> <maxT> <iSave>
```

#### Parameters:

- `<lammps_input_file>`: Path to the LAMMPS input script.
- `<maxT>`: Total number of simulation timesteps.
- `<iSave>`: Interval for saving and publishing simulation data.

#### Example:
```bash
mpirun -np 4 ./singleCell in.lmp4cell 10000 100
```
This command runs the simulation using 4 MPI processes, executes 10,000 timesteps, and saves data every 100 timesteps.