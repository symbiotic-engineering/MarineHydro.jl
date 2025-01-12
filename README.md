![Test Status](https://github.com/symbiotic-engineering/MarineHydro.jl/actions/workflows/run_tests.yml/badge.svg)
![GitHub](https://img.shields.io/github/license/symbiotic-engineering/MDOcean)

### Fully Differentiable Boundary Element Solver for Hydrodynamic Sensitivity Analysis of Wave-Structure Interactions

#### Authors

**Kapil Khanal<sup>a,b</sup>, Carlos A. Michelén Ströfer<sup>b</sup>, Matthieu Ancellin<sup>c</sup>, Maha Haji<sup>a</sup>**

#### Affiliations
- **<sup>a</sup>Cornell University**  
  Ithaca, NY 14850, USA
- **<sup>b</sup>Sandia National Laboratories**  
  Albuquerque, NM 87123, USA
- **<sup>c</sup>Eurobios Mews Labs**  
  Paris, France

## 🌟 Research Highlights

- 📚 **Derivation and discussion** of the discrete adjoint method for the boundary integral equations.
- 💻 **Review and implementation** of a differentiable boundary element solver for marine hydrodynamics in Julia.
- 🌊 **Exact gradient calculation** for a pair of floating hemispheres with respect to their dimensions, separation distance, and wave environment.
- ⚡ **Mechanical power optimization** for a pair of wave energy converters using exact gradients.



Fully-differentiable boundary element solver for marine hydrodynamics. This new solver implements both direct and indirect boundary element formulations and uses two green's function expression, Wu et al, and Delhommeau's varying in their accuracy and speed. 
<img width="632" alt="image" src="https://github.com/user-attachments/assets/16247838-770b-480d-9f2f-d4b0a02054bf" />


> ⚠️ **Note**: This package is **work in progress** 🚧 and a separate public release of the package will be done in the future. This current state of the package contains necessary code to replicate the paper 📄.
>

### Repository: MarineHydro.jl

#### 📂 Folder Structure

- **📁 .github/workflows**  
  Contains workflow files for automated tasks, such as continuous integration (CI).

- **📊 paper**  
  Includes plots and data generated for the paper.

- **📜 src**  
  Source code files for the `BEM.jl` package, including the main functionality.

- **🧪 test**  
  Contains test files and resources to verify the functionality of the source code.

---

### 🚀 How to Run the Code

1. **Install Julia**  
   Ensure you have Julia installed on your system. You can download it from the [JuliaLang website](https://julialang.org/downloads/).

2. **Clone the Repository**  
   Open a terminal and run:  
   ```bash
   git clone https://github.com/symbiotic-engineering/MarineHydro.jl.git
   cd MarineHydro.jl


3. **Install Dependencies**  
   Start Julia from the terminal in the project directory and run the following:  
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
4. **Using the BEM Module**  
Load the module in your Julia session:  
```julia
using BEM
using PyCall
# import your capytaine mesh
cpt = pyimport("capytaine")
radius = 1.0 #fixed
resolution = (10, 10)
cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=resolution) 
cptmesh.keep_immersed_part(inplace=true)

# declare it Julia mesh
mesh = Mesh(cptmesh)
ω = 1.03
ζ = [0,0,1] # HEAVE: will be more verbose in future iteration. define it again even if defined in Capytaine.
F = DiffractionForce(mesh,ω,ζ)
A,B = calculate_radiation_forces(mesh,ζ,ω)
