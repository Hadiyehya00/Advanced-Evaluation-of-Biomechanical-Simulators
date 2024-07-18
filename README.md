# CGoGN

CGoGN is a powerful C++ library designed for mesh data structures, geometric modeling, and geometry processing algorithms. It offers versatility in data representation and processing capabilities, essential for applications in computer graphics and computational geometry.

---

## Features

- **Mesh Data Structures:**
  - Combinatorial maps (1D, 2D, 3D)
  - Graphs (1D)
  - Incidence graphs (1D, 2D)

- **Geometry Processing:**
  - Surface remeshing, subdivision, decimation
  - Surface geometry filtering, curvature computation
  - Surface skeleton extraction
  - Surface deformation
  - Convex hull
  - Hex mesh generation

- **User Interface:**
  - *ImGui* interface with modules and signals
  - *OpenGL* rendering for points, lines, surface, volume meshes

---

## Building CGoGN

To build CGoGN, follow these steps:

1. Create a build directory and navigate into it.
2. Run `cmake <cgogn_path>` to configure the build.
3. Optionally, specify build options such as `-DCMAKE_BUILD_TYPE=""`.
4. On **Linux** and **MacOS**, execute `make -jN` (where N is the number of threads).
5. On **Windows**, use **Visual Studio (2013 or later)** for compilation.

---

## Contribution Guidelines

If you wish to contribute to CGoGN, follow these steps:

1. Fork the GitHub repository.
2. Clone your forked repository to your local machine.
3. Checkout the "develop" branch.
4. Create a new branch named after your username to commit your changes.
5. Push your branch and commits to your GitHub repository.
6. Create a Pull Request from your branch to the "develop" branch of the main repository.

---

## Professional Development in Surgical Simulation

### Problem Statement

Current surgical training methods have many limitations. For example, existing simulators, whether physical mannequins or basic simulation software, are not realistic enough or capable of reproducing all the situations that surgeons may face in an operating room. Moreover, most simulators generally do not allow direct real-time interaction with the mannequin or the handling of complications. This unrealistic and non-interactive training poses an increased risk to the patient due to the inadequate preparation of surgeons.

### Solution Approach

InSimo, a company based in Strasbourg and associated with Inria, designs simulation software dedicated to medical learning through surgery. The company aims to give medical school students their first experiences in a hospital environment, following the principle: "Never the first time on the living."

![Surgical Simulation](https://github.com/Hadiyehya00/Advanced-Evaluation-of-Biomechanical-Simulators/blob/main/cgogn.png)


