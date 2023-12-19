# LatTOp
## Multifuncional Topology Optimisation code, written in Julia, for the optimisation of generic as well as lattice structures.  Implemented goal functions are: Thermal Conductivity, Stiffness, and Phase Change compliance minimisations.
#### This work was initiated by Stefano Piacquadio (stefano.piacquadio@rwth-aachen.de) and developed by Samuel Hayden (samohayden@gmail.com) in the framework of Stefano Piacquadio's research activities at the Institute for Structural Mechanics and Lightweight Design of the RWTH Aachen University

The code performs single and multi-functional topology optimisations of, so far 2D domains. The SIMP method is used. To perform a multifunctional optimisation a Pareto-optimal structure is sought with help of a bisection method.

### Use
1. Generate the domain, mesh and boundary conditions. The code takes in input .msh files. These can be generated with help of the open-source software GMSH https://gmsh.info/.
2. Open the main.jl file and change the path of the mesh file (variable meshpath)
3. Change the filter radius value, depending on the mesh size chosen.
4. Run the main.jl file.
5. Run the optimiseFunctions.jl file, in order to initialise the functions.
6. Choose if you want to use lattice structures, or if ageneric TO should be performed. You can do so by input of the variable `useLattice=true` for lattices, or false.
7. If you want to use lattices, you must define the unit cell topology. For this enter the variable `celltype="bcc"` or `celltype="f2ccz"`. So far only these unit cells were implemented.
8. The single functional optimisation can be started by calling the respective functions, namely `optimiseTopologyThermal(volFrac,"savename",true)` or `optimiseTopologySolid(volFrac,"savename",true)`. The variable volFrac defines the domain's relative density. A savename can be entered and if true, the output will be saved in a .vtk Paraview file.
9. If the Pareto-optimal is sought, the Utopia point must be first found. For this, run the single functional optimisation, separately and save the results in a variable, i.e. `therm_opt=optimiseTopologyThermal(DomainVolumeFraction,"thermal",false)`;`solid_opt=optimiseTopologySolid(DomainVolumeFraction,"solid",false)`. Then, find the optimal weight via: `pareto_opt=getParetoOptimalStructure((therm_opt[2],solid_opt[2]),DomainVolumeFraction)`.
10. The optimal geometry is found by running the multifunctional topology optimisation and setting as a weight the optimal weight found from the procedure above: `optimiseTopology(DomainVolumeFraction,pareto_opt[4],"savename",true)`
11. The result can be opened in Paraview.

