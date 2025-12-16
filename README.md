
# Buoyancy-modified Wilcox (2006) k–omega Model for Incompressible Flows

## Overview

This repository provides an implementation of the **Wilcox (2006)k–omega turbulence model** for incompressible flows with buoyancy, using the **Boussinesq approximation**.
The formulation follows Wilcox (2006) exactly, with buoyancy effects added consistently to the momentum and turbulence transport equations.

The implementation is intended for fundamental studies of buoyancy-driven turbulence and for benchmark flows relevant to natural and mixed convection.

---

The current code is built on OpenFOAM-v2412. 
Below is a brief guide on how to use it.

## 1. Install OpenFOAM 
Follow the guidelines provided at [https://www.openfoam.com/] to install OpenFOAM. 

Note: This code is not compatible with the OpenFOAM Foundation release (often labeled as v8, ..., v12) available at [openfoam.org]

## 2. Set Up Your User Directory
After installation, update your personal "~/.bashrc" file to include OpenFOAM environment variables refer to the installation guidelines. 

e.g. source /usr/lib/openfoam/openfoam2412/etc/bashrc

Then, create the necessary directories for running OpenFOAM by executing:

  > mkdir -p $WM_PROJECT_USER_DIR/{run,applications,src}

Next, overwrite the provided source code (folders named "run" and "applications") by copying the provided files into the appropriate locations.

## 3. Compile the Solver
In a Linux environment, compile the solver located in "applications/mybuoyantBoussinesqSimpleFoam" by running the following commands sequentially:

  > run

  > cd ../applications/buoyantBoussinesqSimpleFoam_kOmega2006
  
  > wclean
  
  > wmake

## 4. Run the Examples  
The "run" folder contains some example problems. 

You can run these examples just like any conventional OpenFOAM tutorial. 

Simply navigate to the example directory and execute:
    
  > sh ./Allclean
  
  > sh ./Allrun
 


---
## Governing equations

This solver considers incompressible flow in which buoyancy is represented by the Boussinesq approximation.  
The Reynolds-averaged momentum and temperature transport equations are

$\frac{\partial U_i}{\partial t} + U_j \frac{\partial U_i}{\partial x_j} = - \frac{\partial P}{\partial x_i} - g_i b (T - T_0) + \frac{\partial}{\partial x_j} \left( \nu \frac{\partial U_i}{\partial x_j} - \overline{u_i u_j} \right)$

$\frac{\partial T}{\partial t} + U_i \frac{\partial T}{\partial x_i} = Q + \frac{\partial}{\partial x_i} \left( a \frac{\partial T}{\partial x_i} - \overline{\theta u_i} \right)$

where $U_i$ and $u_i$ are the mean and fluctuating velocity components, $T$ and $\theta$ are the mean and fluctuating temperature, $x_i$ denotes the Cartesian coordinates, and $t$ is time.  
$P$ is the mean kinematic pressure, $g_i$ is the gravitational acceleration, $b$ is the thermal expansion coefficient, $T_0$ is a reference temperature, $\nu$ is the kinematic viscosity, $Q$ is the volumetric heat source per unit heat capacity, and $a$ is the thermal diffusivity.

When a lower wall is present, the value of $T_0$ merely shifts the hydrostatic reference state and does not affect the fluid motion; for simplicity, $T_0$ is set to zero.  
The overbar denotes Reynolds averaging, and $\overline{u_i u_j}$ and $\overline{\theta u_i}$ are the Reynolds stress tensor and turbulent heat flux vector, respectively.

---

## Turbulence closure

The Reynolds stress tensor and turbulent heat flux are modelled using the eddy-viscosity and gradient-diffusion hypotheses:

$\overline{u_i u_j} = \frac{2}{3} k \delta_{ij} - 2 \nu_T S_{ij}, \qquad \overline{\theta u_i} = - a_T \frac{\partial T}{\partial x_i}$

with

$S_{ij} = \frac{1}{2} \left( \frac{\partial U_i}{\partial x_j} + \frac{\partial U_j}{\partial x_i} \right), \qquad a_T = \frac{\nu_T}{\mathit{Pr}_T}$

Here $S_{ij}$ is the mean strain-rate tensor, $\nu_T$ is the turbulent viscosity, $a_T$ is the turbulent thermal diffusivity, and the turbulent Prandtl number is set to $\mathit{Pr}_T = 0.89$.

---

## Wilcox (2006) k–omega model with buoyancy

The notation for the model variables and constants follows the book of Wilcox (2006), with only buoyant terms added.

The turbulent viscosity is defined as

$\nu_T = \frac{k}{\tilde{\omega}}, \qquad \tilde{\omega} = \max \left\lbrace \omega, C_{lim} \sqrt{ \frac{2 S_{ij} S_{ij}}{\beta^{\ast}} } \right\rbrace$

The transport equations for $k$ and $\omega$ are

$\frac{\partial k}{\partial t} + U_j \frac{\partial k}{\partial x_j} = \mathcal{P} + \mathcal{P}_b - \beta^{\ast} k \omega + \frac{\partial}{\partial x_j} \left[ \left( \nu + \sigma^{\ast} \frac{k}{\omega} \right) \frac{\partial k}{\partial x_j} \right]$

$\frac{\partial \omega}{\partial t} + U_j \frac{\partial \omega}{\partial x_j}= \alpha \frac{\omega}{k} \left( P + C_{\omega b}^{+} \max( P_b,0)  + C_{\omega b}^{-} \min( P_b,0)  \right)- \beta_0 f_{\beta} \omega^2 + \frac{\sigma_d}{\omega} \frac{\partial k}{\partial x_j} \frac{\partial \omega}{\partial x_j} + \frac{\partial}{\partial x_j} [ ( \nu + \sigma \frac{k}{\omega} ) \frac{\partial \omega}{\partial x_j} ]$

The production terms due to the mean velocity gradient and buoyancy are

$\mathcal{P} = - \overline{u_i u_j} \frac{\partial U_i}{\partial x_j}, \qquad \mathcal{P}_b = - g_i b \, \overline{\theta u_i}$

---

## Auxiliary relations

$\sigma_d = 0$   for   $\frac{\partial k}{\partial x_j} \frac{\partial \omega}{\partial x_j} \le 0$,
$\sigma_d = \sigma_{do}$ for $\frac{\partial k}{\partial x_j} \frac{\partial \omega}{\partial x_j} > 0$

$f_{\beta} = \frac{1 + 85 \chi_\omega}{1 + 100 \chi_\omega}, \qquad \chi_\omega = \left| \frac{ \Omega_{ij} \Omega_{jk} S_{ki} }{ (\beta^{\ast} \omega)^3 } \right|$

with

$\Omega_{ij} = \frac{1}{2} \left( \frac{\partial U_i}{\partial x_j} - \frac{\partial U_j}{\partial x_i} \right)$

---

## Model constants

The closure coefficients are $C_{lim} = 0.875$, $\beta^{\ast} = 0.09$, $\sigma^{\ast} = 0.6$, $\alpha = 0.52$, $\beta_0 = 0.0708$, $\sigma = 0.5$, and $\sigma_{do} = 0.125$.  
For buoyancy effects, the default values $C_{\omega b}^{+} = 1$ and $C_{\omega b}^{-} = -2$ are adopted, consistent with the implementation in ANSYS CFX.

---

## Wall boundary conditions

At no-slip walls, the boundary conditions are

$k = 0, \qquad \omega = \frac{6 \nu}{\beta_0 n^2}$

where, for a cell-centred discretisation of $\omega$, $n$ is the wall-normal distance from the wall to the centre of the first off-wall grid cell.

For user convenience, if the boundary condition for $\omega$ is specified as "type fixedValue", regardless of the prescribed value, the code automatically replaces it with this boundary condition and applies it accordingly.


---


## Example Cases

The repository includes example cases for:

* Rayleigh–Bénard convection
* Vertical heated convection
* Unstably stratified Couette flow
* Stably stratified plane channel flow
* Internally heated convection

see https://arxiv.org/abs/2512.01308

---
