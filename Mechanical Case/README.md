# Mechanical Case - Observer Design

This folder contains the implementation of observer design for mechanical port-Hamiltonian systems.

In order to run the Simulink scheme before run the code: `Init_Simulink_param.m` .

The Simulink scheme presented inside the folder `nonlinear_mech_observer` have an observer gain that can be changed according to the constant one obtained form the output of `MechCode.m` or other values is order to appreciate the different results on the other side no scheduled gain is implemented in those Simulink scheme.

`MechCode.m` is the code that implement the constant and scheduled observer gain procedure, the code produces also the same images used in the paper.

## Contents

- MATLAB/Simulink simulation files
- Observer implementation code
- Test cases and validation results

## Description

Implementation of gain-scheduled observer using LMI-based synthesis for mechanical port-Hamiltonian systems, including polytopic LPV embedding and Lyapunov stability analysis.
