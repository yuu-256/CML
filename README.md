# Coupled Map Lattice (CML) model for cloud dynamics

### Usage
1. Make random initial conditions for the CML model by running `utils/mkini.py`.  
2. Set configuration file path as initial conditions by editing `code/main.f90`.  
3. Set the number of steps and simulation parameters in `code/main.f90`.
3. Run `./run.sh` and it will compile and run the code, and plot the results(edit utils/plot_2D.py for different plots).  

### Structure  
code/main.f90: main program  

### Debugging  
NaN happens about 1000 steps after the start. 
x = 47, 48  

**checked:**  
subroutine cml_sim 
subroutine buoyancy_dragging  
subroutine viscosity_pressure 

**not yet:** 
subroutine diffusion 
subroutine phase_transition 
subroutine lagrangian 

Implementations of each differential equation seem to have problems (results seem to diverge?). 
Conservation laws are not satisfied
