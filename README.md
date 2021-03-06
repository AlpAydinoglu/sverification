Here, you can find the code for the paper 'Stability Analysis of Complementarity Systems with Neural Network Controllers' from the conference HSCC 2021.

# Dependencies

In order to run the examples, both MATLAB and Python are required.

## MATLAB
The linear complementarity problems (LCPs) generated by this library are solved using [http://pages.cs.wisc.edu/~ferris/path.html](PATH). 

MOSEK is used in order to solve the optimization problems (https://www.mosek.com/).

Optimization problems are formulated using https://yalmip.github.io/ (YALMIP).

`pathlcp`, `mosek`, `yalmip` will need to be in the MATLAB path for the examples to run.

## Python

Python3 is used for all of the examples with the following additional packages:

`numpy`: https://numpy.org/

`pytorch`: https://pytorch.org/

`scipy.io`: https://docs.scipy.org/doc/scipy/reference/tutorial/io.html

`numpy`, `pytorch`, `scipy.io` should be imported for the examples to run.

## Functionality

The library can be used to verify stability of linear complementarity systems with neural network controllers.

## Examples

`double_int`: Double integrator (Example 5.1)

`cartpole`: Cart-pole with soft walls (Example 5.2)

`box_friction`: Box with friction (Example 5.3)

`five_carts`: Five carts with soft contacts (Example 5.4)
