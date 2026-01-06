# Numerical Methods

This document describes the **numerical strategy actually used in the code**, including integrators, termination criteria, and practical heuristics.  

---

## 1. Phase-space formulation

All curves—both null rays and timelike observer worldlines—are evolved in **phase space**

$$
(x^\mu, p_\mu),
$$

using Hamilton’s equations derived from

$$
H(x,p) = \tfrac12\, g^{\mu\nu}(x)\, p_\mu p_\nu.
$$

This formulation has several practical advantages:

- A single evolution system applies to both null and timelike curves.
- The equations remain regular at the future horizon in ingoing coordinates.
- The velocity update $dx^\mu/d\lambda = g^{\mu\nu} p_\nu$ is explicit and simple.

---

## 2. Metric derivatives

The contravariant metric $g^{\mu\nu}(r,\theta)$ is implemented analytically.

Its derivatives with respect to coordinates are handled as follows:

- $\partial_r g^{\mu\nu}$ and $\partial_\theta g^{\mu\nu}$ are computed numerically.
- $\partial_v g^{\mu\nu} = 0$ and $\partial_\phi g^{\mu\nu} = 0$, reflecting stationarity and axisymmetry.

These derivatives enter only through the momentum evolution equation

$$
\frac{dp_\mu}{d\lambda}
= -\tfrac12 \, \partial_\mu g^{\alpha\beta} \, p_\alpha p_\beta.
$$

---

## 3. Time stepping and integrator

All geodesics are evolved using a **fixed-step classical fourth-order Runge–Kutta (RK4)** scheme.

- No adaptive error control is employed.
- Stability is achieved through conservative step sizes and geometric regularity of the coordinate system.

### Step-size heuristic

For null rays, the affine step size $\Delta\lambda$ is scaled heuristically with radius:
- Smaller steps are used near the black hole,
- Larger steps are used far from the hole.

This heuristic is purely practical and is not derived from an error estimator.

For timelike trajectories, integration is performed in proper time $\tau$ with a user-specified step size.

---

## 4. Timelike geodesic integration

Timelike observer trajectories are evolved using the same Hamiltonian system as null rays, with additional normalization steps:

1. The initial four-velocity $u^\mu$ is normalized so that

   $$
   g_{\mu\nu} u^\mu u^\nu = -1.
   $$

2. The corresponding momentum is set as

   $$
   p_\mu = g_{\mu\nu} u^\nu.
   $$

3. After each RK4 step:
   - $u^\mu$ is recomputed via $u^\mu = g^{\mu\nu} p_\nu$,
   - $u^\mu$ is renormalized to control numerical drift,
   - $p_\mu$ is reconstructed from the normalized velocity.

This procedure sacrifices exact symplecticity in favor of numerical robustness for long timelike integrations.

---

## 5. Camera tetrad construction

At each observer event, an orthonormal tetrad is constructed numerically:

- The timelike leg is fixed to the normalized observer four-velocity.
- The spatial legs are constructed via Gram–Schmidt orthonormalization with respect to the Kerr metric.

A sign check ensures that the designated “forward” spatial leg points toward decreasing radial coordinate $r$.

Yaw, pitch, and roll are implemented as explicit rotations of the spatial triad prior to orthonormalization.

The tetrad defines the mapping between pixel directions and spacetime directions.

---

## 6. Null ray initialization

For each pixel:

- A normalized direction vector is constructed in the camera frame.
- This direction is combined with the timelike leg to form a **past-directed null vector**.
- The null vector is mapped to spacetime coordinates using the tetrad.
- The covariant momentum $p_\mu = g_{\mu\nu} k^\nu$ is used as the initial condition.

Absolute scaling of the null momentum is irrelevant and is not fixed.

---

## 7. Termination conditions

Null rays are evolved until one of the following occurs:

### Escape
If the radial coordinate exceeds a prescribed threshold $r_{\text{escape}}$, the ray is considered to have reached the asymptotic region.

The angular direction at escape is recorded and used to sample the sky texture.

### Horizon crossing
If the ray crosses the future event horizon $r_+$, it is classified as captured.

Crossing is detected by a sign change in $r - r_+$ between integration steps.

---

## 8. Near-horizon filtering

For cameras outside the horizon, additional logic suppresses numerical artifacts near $r_+$:

- If a ray’s minimum radius during integration falls below a small tolerance above the horizon,
  it is forcibly classified as captured.

This filter is a pragmatic measure to reduce visual speckling caused by finite step sizes.

---

## 9. Inside-horizon classification

When the camera is inside the horizon, rays cannot escape in the usual sense.

In this regime, rays are classified based on whether they ever exhibit outward radial motion during integration.

This heuristic produces visually consistent images but should not be interpreted as a physically invariant classification.

---

## 10. Angular singularities and coordinate issues

The coordinate system contains the usual Kerr singularities:

- The axis $\theta = 0, \pi$ introduces $\sin\theta$ factors in metric components.
- No special regularization is applied at the axis.

Rendering setups are expected to avoid exact axis alignment.

---

## 11. Sky sampling and color handling

The sky is represented by a static equirectangular texture.

- $\phi$ maps to the horizontal texture coordinate,
- $\theta$ maps to the vertical coordinate.

No intensity modulation, redshift, or Doppler factor is applied.

A simple gamma correction is applied to the sampled RGB values.

---

## 12. Limitations and scope

The numerical scheme prioritizes **clarity and robustness** over optimal performance or strict conservation.

Known limitations include:
- lack of adaptive step-size control,
- absence of conserved-quantity monitoring,
- heuristic near-horizon filtering,
- no radiative transfer or frequency effects.

These choices reflect the current scope of the project and may be revisited in future extensions.

---

## 13. Intended use

The numerical methods implemented here provide a **clean, geometry-aware framework** for relativistic ray tracing in Kerr spacetime.

They are suitable for qualitative and semi-quantitative studies of strong-field lensing and observer-dependent imaging, but are not intended as high-precision geodesic solvers.
