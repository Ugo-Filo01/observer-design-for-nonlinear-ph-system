# Electromechanical Case - Observer Design

This folder contains the implementation of observer design for electromechanical port-Hamiltonian systems.

The script `ElectromechCode.m` defines both the constant-gain and gain-scheduled procedures, following the methodology presented in the thesis.

The folder `linearized_electromech_observer` contains the Simulink model associated with the linearized version of the observer.  
Before running the Simulink simulation, execute the script `Init_Simulink_param_lin_Electromech.m`.  
The observer gain can be adjusted directly within the Simulink model.

**Additional note:** the linearization of the mechanical model is not included, as the approach follows the framework proposed in the paper by Spirito, where only the electromechanical case fits one of the two classes presented.

The folder `nonlinear_electromech_observer` contains the Simulink model corresponding to the nonlinear implementation of the observer.  
Before running the Simulink simulation, execute the script `Init_Simulink_param_nl_Electromech.m`.
## Contents

- MATLAB/Simulink simulation files

- Observer implementation code

- Test cases and validation results
