# Conventions

This document records the **coordinate, sign, and normalization conventions used throughout the code**.  
It exists solely to eliminate ambiguity; it contains no derivations.

---

## Coordinates and ordering

All spacetime quantities are expressed in **ingoing (advanced) Kerr coordinates**

$$
x^\mu = (v,\ r,\ \theta,\ \phi),
$$

with this ordering used **everywhere** in the code and Python bindings.

- $v$: advanced time coordinate  
- $r$: Boyer–Lindquist–like radial coordinate  
- $\theta$: polar angle  
- $\phi$: azimuthal angle  

Arrays representing spacetime points or vectors always follow this ordering.

---

## Metric signature and units

- Metric signature: $(-, +, +, +)$
- Units: $G = c = 1$

All quantities are dimensionless or expressed in geometric units set by the mass parameter $M$.

---

## Index conventions

- Greek indices $\mu,\nu,\ldots$ run over spacetime coordinates $(v,r,\theta,\phi)$.
- Parenthesized indices $(a),(b),\ldots$ label **tetrad-frame components**.

Summation over repeated indices is assumed.

---

## Covariant vs contravariant quantities

- Positions are represented as $x^\mu$.
- Momenta are represented as **covariant** components $p_\mu$.
- Velocities and tetrad legs are typically contravariant vectors.

Raising and lowering of indices is performed using the Kerr metric:

$$
p_\mu = g_{\mu\nu} u^\nu, \qquad
u^\mu = g^{\mu\nu} p_\nu.
$$

---

## Geodesic parameterization

- **Null rays** are evolved with respect to an affine parameter $\lambda$.
- **Timelike observer trajectories** are evolved with respect to proper time $\tau$.

The same Hamiltonian system is used in both cases; the curve type is determined by the initial normalization.

---

## Normalization conditions

- Timelike four-velocities satisfy

  $$
  g_{\mu\nu} u^\mu u^\nu = -1.
  $$

- Null momenta satisfy

  $$
  g^{\mu\nu} p_\mu p_\nu = 0.
  $$

Timelike trajectories are explicitly renormalized during integration to control numerical drift.

---

## Tetrad conventions

- $e_{(0)}^\mu$: timelike tetrad leg, aligned with the observer four-velocity.
- $e_{(1)}^\mu$: camera “forward” direction.
- $e_{(2)}^\mu$: vertical image-plane direction.
- $e_{(3)}^\mu$: horizontal image-plane direction.

The spatial triad is orthonormalized with respect to the Kerr metric.

A sign check enforces that the forward direction corresponds to **decreasing radial coordinate $r$**.

---

## Ray-tracing direction

Null rays are **back-traced** from the camera into the external sky:

- Rays are launched from the camera event.
- Integration proceeds outward in spacetime until escape or horizon crossing.
- The absolute scaling of the null momentum is irrelevant and not fixed.

---

## Horizon definition

The future event horizon is defined as

$$
r_+ = M + \sqrt{M^2 - a^2}.
$$

Horizon crossing is detected by a sign change in $r - r_+$ between integration steps.

---

## Angular conventions for sky mapping

- $\phi \in [0, 2\pi)$ maps to the horizontal sky texture coordinate.
- $\theta \in [0, \pi)$ maps to the vertical sky texture coordinate.

The sky texture is assumed to be equirectangular.

---

## Summary

All conventions listed here are **fixed by the implementation** and should be treated as part of the model definition.  
Any extensions or refactors should update this document accordingly.
