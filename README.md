# ODESolver
Practicing C++ by implementing various ODE Solvers. Used in other projects.

## Context 
To solve differential equations/initial value problems, an ODE solver must be used (e.g., `ode45` in Matlab and `solve_ivp` in Python/scipy). I wanted to practice my C++ and write my own ODE solver which will be used in future projects modeling the dynamics of systems over time. 


I took an OOP approach to allow for the easy implementation of new solvers. Currently, an abstract ODE Solver is used (so that other systems can use its type and reference the methods) which is derived with two classes: the `ForwardEuler` and `RKF45` class. The Forward Euler is a naive approach to integration and simply adds the previous state to the product of the time-difference and state derivative. The error from this method will blow up quickly.
RK-45, on the other hand, is the workhorse of integration schemes and is able to compute integrations with 5th order local error and 4th order global error with relatively little work. This is used in `ode45` and the default method in `solve_ivp`.

## Includes
Requires the `Eigen` C++ library which can be found [here]().
Also uses the `MatplotLib` C++ header to plot data directly in C++ which can be found [here]().
`CMake` and `make` are used also used to build the program.

## Running 
To run, clone the github in the directory of your choosing via 
```
git clone https://github.com/gagandeepthapar/ODESolver/tree/main
```


Ensure that the required libraries are installed. My system adds the `Eigen` and `MatplotLib` headers directly in my include path so I don't have to explicity call it. You can either add the libaries in the same project directory or where you'd like (e.g., `usr/local/include`) and add them to your include path:
```
$ export CPLUS_INCLUDE_PATH="usr/local/include:$CPLUS_INCLUDE_PATH"
```


Once installed, you should be able to run the following commands in the project directory to setup the build directory and output, respectively:
```
$ cmake -S . -B build
$ make -C build
```


The output file should be called `example.out` which can be renamed directly in the `CMakeLists.txt` file.
