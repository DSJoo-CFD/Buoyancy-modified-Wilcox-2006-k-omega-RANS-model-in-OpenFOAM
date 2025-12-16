
# Buoyancy-modified Wilcox (2006) (k)–(\omega) Model for Incompressible Flows

## Overview

This repository provides an implementation of the **Wilcox (2006) (k)–(\omega) turbulence model** for incompressible flows with buoyancy, using the **Boussinesq approximation**.
The formulation follows Wilcox (2006) exactly, with buoyancy effects added consistently to the momentum and turbulence transport equations.

The implementation is intended for fundamental studies of buoyancy-driven turbulence and for benchmark flows relevant to natural and mixed convection.

---

## Governing Equations

### Reynolds-averaged equations

Buoyancy is modelled using the Boussinesq approximation. The Reynolds-averaged momentum and temperature equations are

[
\frac{\partial U_i}{\partial t}

* U_j \frac{\partial U_i}{\partial x_j}
  = - \frac{\partial P}{\partial x_i}

- g_i b (T - T_0)

* \frac{\partial}{\partial x_j}
  \left( \nu \frac{\partial U_i}{\partial x_j} - \overline{u_i u_j} \right),
  ]

[
\frac{\partial T}{\partial t}

* U_i \frac{\partial T}{\partial x_i}
  = Q
* \frac{\partial}{\partial x_i}
  \left( a \frac{\partial T}{\partial x_i} - \overline{\theta u_i} \right).
  ]

Here, (U_i) and (T) denote the mean velocity and temperature, respectively.
The reference temperature (T_0) is set to zero, as it only shifts the hydrostatic pressure in wall-bounded flows.

---

## Turbulence Closure

### Reynolds stress and turbulent heat flux

The Reynolds stress tensor and turbulent heat flux are modelled using the eddy-viscosity and gradient-diffusion hypotheses:

[
\overline{u_i u_j} = \frac{2}{3}k\delta_{ij} - 2\nu_T S_{ij},
\qquad
\overline{\theta u_i} = - a_T \frac{\partial T}{\partial x_i},
]

with
[
S_{ij} = \frac{1}{2}
\left(
\frac{\partial U_i}{\partial x_j} +
\frac{\partial U_j}{\partial x_i}
\right),
\qquad
a_T = \frac{\nu_T}{\mathit{Pr}_T}.
]

The turbulent Prandtl number is fixed at
[
\mathit{Pr}_T = 0.89.
]

---

## (k)–(\omega) Model (Wilcox 2006)

The turbulent viscosity is defined as
[
\nu_T = \frac{k}{\tilde{\omega}}, \qquad
\tilde{\omega} = \max \left(
\omega,,
C_{lim} \sqrt{\frac{2S_{ij}S_{ij}}{\beta^{\ast}}}
\right).
]

The transport equations for turbulent kinetic energy (k) and specific dissipation rate (\omega) are

[
\frac{\partial k}{\partial t}

* U_j \frac{\partial k}{\partial x_j}
  =
  \mathcal{P} + \mathcal{P}_b
  -\beta^{\ast} k \omega
* \frac{\partial }{\partial x_j}
  \left[
  \left( \nu + \sigma^{\ast} \frac{k}{\omega} \right)
  \frac{\partial k}{\partial x_j}
  \right],
  ]

[
\frac{\partial \omega}{\partial t}

* U_j \frac{\partial \omega}{\partial x_j}
  =
  \alpha \frac{\omega}{k}
  \left( \mathcal{P} + C_{\omega b} \mathcal{P}*b \right)
  -\beta_0 f*{\beta} \omega^{2}
* \frac{\sigma_d}{\omega}
  \frac{\partial k}{\partial x_j}
  \frac{\partial \omega}{\partial x_j}
* \frac{\partial }{\partial x_j}
  \left[
  \left( \nu + \sigma \frac{k}{\omega} \right)
  \frac{\partial \omega}{\partial x_j}
  \right].
  ]

The production terms are
[
\mathcal{P} = -\overline{u_i u_j}
\frac{\partial U_i}{\partial x_j},
\qquad
\mathcal{P}_b = - g_i b \overline{\theta u_i}.
]

The buoyancy contribution in the (\omega)-equation is controlled by the model constant (C_{\omega b}).

---

## Model Constants

All closure coefficients follow Wilcox (2006):

* (C_{lim} = 0.875)
* (\beta^{\ast} = 0.09)
* (\sigma^{\ast} = 0.6)
* (\alpha = 0.52)
* (\beta_0 = 0.0708)
* (\sigma = 0.5)
* (\sigma_{do} = 0.125)

---

## Wall Boundary Conditions

At no-slip walls, the boundary conditions are
[
k = 0,
\qquad
\omega = \frac{6\nu}{\beta_0 n^2},
]
where (n) is the wall-normal distance to the centre of the first off-wall cell.

---

## Example Cases

The repository includes example cases for:

* Rayleigh–Bénard convection
* Vertical heated convection
* Unstably stratified Couette flow
* Stably stratified plane channel flow
* Internally heated convection

---
* **OpenFOAM solver 구조 설명 섹션**
* **논문 Appendix ↔ README 대응 표**
