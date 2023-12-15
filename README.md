# Thesis
Reposity for all Thesis related code

# test git
 
# File Structure
- src: Contains the main scripts I am working on

- res: Contains sub-folders with resources I will use
- res/meshes: Contains the meshes and .geo files I have created and can use.
- res/results: This is gitignored, but this is where results will be stored.

- misc: Contains subfolders of misc data.
- misc/Old Code: Contains old functional code I have written that I may need to reference.
- misc/Reference Code: Contains old code using libraroes that I may need to reference.

# Running the Script
1. Execute "src/main.jl" in Julia REPL
2. run commands in REPL terminal
	- As of writting this, there are multiple "solve***()" functions which take in an arguement. The argument does nothing, but it will be used as the porosity variable (and potentially others) once I impliment the material properties of the cells.
	- The optimiseTopology function optimises the thermal problem. This still needs to be generalised#LatTOp
