# Cascade and Parallel PID Tuning for Inverted Pendulum

This repository contains MATLAB codes for tuning and simulating PID controllers for an inverted pendulum on a cart. The project focuses on two control structures: cascade PID control and parallel PID control. The goal is to stabilize the pendulum angle and control the cart position while satisfying performance constraints such as settling time, maximum pendulum angle, and actuator saturation.

## Project Overview

The inverted pendulum is an inherently unstable nonlinear system. In this project, both linear and nonlinear models are considered, and PID-based controllers are designed to stabilize the system. Since manual tuning of cascade and parallel PID gains is time-consuming and highly dependent on trial and error, the provided codes automate the search process and evaluate candidate gains based on predefined performance criteria.

## Control Structures

### Cascade PID Control

In the cascade structure, the outer loop controls the cart position and generates a reference angle for the pendulum. The inner loop then controls the pendulum angle and produces the final control force applied to the system.

### Parallel PID Control

In the parallel structure, two PID controllers work simultaneously. One controller acts on the cart position error, while the other acts on the pendulum angle error. Their outputs are combined to generate the final control input.

## Main Features

- MATLAB implementation of cascade PID tuning
- MATLAB implementation of parallel PID tuning
- Fixed-step ode4 simulation with 0.001 s step size
- Linear and nonlinear inverted pendulum models
- Force saturation handling
- Settling time and overshoot evaluation
- Maximum pendulum angle checking
- Automatic search for acceptable PID gains
- Plotting of cart position, pendulum angle, controller output, and applied force

## Performance Criteria

The tuning process evaluates controllers based on the following conditions:

- Cart position settling time
- Pendulum angle settling time
- Maximum pendulum angle during simulation
- Maximum applied control force
- Final steady-state error
- Stability of the closed-loop response

## Repository Files

Suggested file structure:

```text
.
├── Cascade_PID_Tuner.m
├── Parallel_PID_Tuner.m
├── README.md
└── results/
```

## How to Run

1. Open MATLAB.
2. Place the MATLAB files in the same working directory.
3. Run the desired tuner file:

```matlab
Cascade_PID_Tuner
```

or:

```matlab
Parallel_PID_Tuner
```

4. The script will search for suitable PID gains and display the final results and plots.

## Notes

The controller gains obtained from the code are intended for the specified model parameters. If the physical parameters of the cart, pendulum, actuator, or saturation limits change, the tuning process should be repeated.

## Author

This project was developed for the simulation and control of an inverted pendulum system using MATLAB and Simulink.
