## ***Advanced-Evaluation-of-Biomechanical-Simulators***

CGoGN is a powerful C++ library designed for mesh data structures, geometric modeling, and geometry processing algorithms. It offers versatility in data representation and processing capabilities, essential for applications in computer graphics and computational geometry.

![CGoGN](https://github.com/Hadiyehya00/Advanced-Evaluation-of-Biomechanical-Simulators/blob/main/cgogn.png)

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
2. Run cmake <cgogn_path> to configure the build.
3. Optionally, specify build options such as -DCMAKE_BUILD_TYPE="".
4. On **Linux** and **MacOS**, execute make -jN (where N is the number of threads).
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

Surgery continues to play a vital role in the treatment of numerous diseases, with its quality and scope improving as progressive development allows not only for saving patients' lives but also greatly enhancing their quality. However, the entirety of surgical practice is precarious and risky. For instance, operating room incidents can arise from undetected anomalies that become apparent or from imperfect techniques. This underscores the need for comprehensive and practical surgical training that prepares young doctors to handle any clinical situation, whether routine or critical. Surgeons should undergo practical training, akin to how pilots are trained on simulators before taking control. Biomechanical simulators that accurately mimic reality represent the future of medical education, ensuring that trainee surgeons can practice both routine and critical cases without endangering patients.

### Solution Approach

![CGoGN](https://github.com/Hadiyehya00/Advanced-Evaluation-of-Biomechanical-Simulators/blob/main/endoscope.png)

*Image used with permission from InSimo.*
*Source: [(https://www.insimo.com/display-bariatric/)]*


Current surgical training methods have numerous shortcomings. For instance, existing simulators, whether physical mannequins or basic simulation software, are neither sufficiently realistic nor capable enough to replicate the full range of situations surgeons may encounter in an operating room. Moreover, most simulators generally do not allow for real-time direct interaction with the mannequin or management of complications. This unrealistic and non-interactive training poses an increased risk to patients due to surgeons' lack of preparedness.

To address these issues, it is crucial to develop an ultra-realistic and interactive biomechanical simulator capable of faithfully recreating the operating room environment. Such a simulator should be able to handle both straightforward and complex scenarios, as well as complications, with real-time direct interaction. Such training is essential to prepare novice surgeons and facilitate their early transition to interacting with real patients.

In this context, the primary challenge to address in this project is the creation of simulation cases for validating solver performance. Specifically, it would be beneficial to compare different test case implementations to determine which offers the best balance between realism, accuracy, and computational time for simulation. This would later enable standardization of testing methods and objective comparison of different implementations, thereby facilitating the selection of the best implementation to enhance surgeons' resolution skills and ensure patient safety.

<div style="text-align: center; margin-bottom: 20px;">
    <img src="https://github.com/Hadiyehya00/Advanced-Evaluation-of-Biomechanical-Simulators/raw/main/insimo.png" alt="InSimo" width="200" style="margin-bottom: 10px;">
    <br>
    <img src="https://github.com/Hadiyehya00/Advanced-Evaluation-of-Biomechanical-Simulators/raw/main/ICube.jpg" alt="ICube" width="160" style="margin-top: 10px;">
</div>

*© 2024 ICube. Tous droits réservés.*

*© 2024 ISimo. Tous droits réservés.*

---
