# Coupled Map Lattice (CML) model for cloud dynamics

Debugging  
NaN happens about 1000 steps after the start. 
x = 47, 48  

checked: 
subroutine cml_sim 
subroutine buoyancy_dragging  
subroutine viscosity_pressure 

not yet: 
subroutine diffusion 
subroutine phase_transition 
subroutine lagrangian 

Implementations of each differential equation seem to have problems (results seem to diverge?). 
Conservation laws are not satisfied
