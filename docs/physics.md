# Physics Model

This document describes the **physical and mathematical model actually implemented in the code**.  
It is intentionally precise and avoids features that are *not* present in the current implementation.

The simulation models **relativistic imaging in Kerr spacetime**, as experienced by a camera moving along a **timelike geodesic**, with images formed by **ray tracing null geodesics** into an external sky.

---

## 1. Spacetime and coordinates

The background spacetime is **Kerr**, characterized by mass $M$ and spin parameter $a$.

All computations are performed in **ingoing (advanced) Eddington–Finkelstein–type coordinates**

$$
x^\mu = (v, r, \theta, \phi),
$$

with metric signature $(-, +, +, +)$.

We define the standard Kerr scalars

$$
\Sigma = r^2 + a^2 \cos^2\theta, \qquad
\Delta = r^2 - 2 M r + a^2.
$$

The use of ingoing coordinates ensures regularity at the **future event horizon**

$$
r_+ = M + \sqrt{M^2 - a^2}.
$$

Both the **covariant metric** $g_{\mu\nu}$ and the **contravariant metric** $g^{\mu\nu}$ are implemented explicitly in the code.

---

## 2. Geodesic equations (Hamiltonian formulation)

All geodesics—both null rays and timelike observer trajectories—are evolved using a **Hamiltonian phase-space formulation**.

The Hamiltonian is

$$
H(x,p) = \tfrac12 \, g^{\mu\nu}(x)\, p_\mu p_\nu,
$$

where $p_\mu$ is the covariant momentum.

Hamilton’s equations are

$$
\frac{dx^\mu}{d\lambda} = \frac{\partial H}{\partial p_\mu}
= g^{\mu\nu}(x)\, p_\nu,
$$

$$
\frac{dp_\mu}{d\lambda} = -\frac{\partial H}{\partial x^\mu}
= -\tfrac12 \, \partial_\mu g^{\alpha\beta}(x)\, p_\alpha p_\beta.
$$

These equations are implemented directly in the C++ core.

### Metric derivatives

The metric depends explicitly only on $r$ and $\theta$.  
Accordingly:
- $\partial_r g^{\mu\nu}$ and $\partial_\theta g^{\mu\nu}$ are computed numerically,
- $\partial_v g^{\mu\nu} = \partial_\phi g^{\mu\nu} = 0$.

---

## 3. Null vs timelike geodesics

The same Hamiltonian system is used for both null and timelike curves.

- **Null rays** satisfy

  $$
  g^{\mu\nu} p_\mu p_\nu = 0.
  $$

  These are used for ray tracing.

- **Timelike observer trajectories** satisfy

  $$
  g^{\mu\nu} p_\mu p_\nu = -1,
  $$

  enforced by explicit normalization during integration.

The character of the geodesic is determined entirely by the **initial momentum**, not by separate equations of motion.

---

## 4. Camera worldline (timelike observer)

The camera follows a **timelike geodesic** $$x^\mu(\tau)$$, with four-velocity

$$
u^\mu = \frac{dx^\mu}{d\tau}.
$$

Implementation details:
- The user supplies an initial position $x^\mu_0$ and four-velocity $u^\mu_0$.
- The four-velocity is normalized so that $g_{\mu\nu} u^\mu u^\nu = -1$.
- The corresponding covariant momentum is

  $$
  p_\mu = g_{\mu\nu} u^\nu.
  $$

- After each integration step, the four-velocity is recomputed from $p_\mu$ and renormalized to control numerical drift.

This procedure yields a stable timelike worldline suitable for frame-by-frame rendering.

---

## 5. Local camera frame (tetrad)

At each camera event, an **orthonormal tetrad**

$$
\{ e_{(0)}^\mu, e_{(1)}^\mu, e_{(2)}^\mu, e_{(3)}^\mu \}
$$

is constructed.

- $e_{(0)}^\mu = u^\mu$ is the timelike leg.
- $e_{(1)}^\mu$ defines the camera “forward” direction.
- $e_{(2)}^\mu$ and $e_{(3)}^\mu$ span the transverse image plane.

The spatial triad is orthonormalized with respect to the Kerr metric.  
A sign check ensures that the forward direction corresponds to **decreasing radial coordinate** (inward-pointing).

Yaw, pitch, and roll are implemented as rotations of the spatial triad prior to orthonormalization.

This tetrad fully defines the camera’s instantaneous rest frame and orientation.

---

## 6. Ray construction and back-tracing

For each image pixel:
1. A direction is constructed in the **camera frame** based on the field of view.
2. This direction is combined with the timelike leg to form a **past-directed null vector** in the tetrad frame.
3. The null vector is mapped to coordinates via

   $$
   k^\mu = e_{(a)}^\mu k^{(a)}.
   $$

4. The corresponding covariant momentum is obtained as

   $$
   p_\mu = g_{\mu\nu} k^\nu.
   $$

The resulting $(x^\mu, p_\mu)$ serves as the initial condition for null geodesic integration.

This procedure implements **ray back-tracing**: rays are launched from the camera and integrated outward to determine what the camera sees.

---

## 7. Ray termination and classification

Null geodesics are integrated until one of the following occurs:

### Escape
If the radial coordinate exceeds a user-defined threshold $r_{\mathrm{escape}}$, the ray is treated as having reached the asymptotic region.

The angular direction $(\theta, \phi)$ at escape is recorded and used to sample the external sky texture.

### Horizon capture
If the ray crosses the future event horizon $r_+$, it is classified as captured and contributes a black pixel.

### Near-horizon filtering
For cameras outside the horizon, additional logic suppresses numerical “speckle” artifacts by marking rays as captured if they approach within a small tolerance of $r_+$.

### Inside-horizon behavior
For cameras inside the horizon, rays are classified based on whether they ever exhibit outward radial motion during integration.

These classification rules are pragmatic and are intended to produce visually consistent images near the horizon.

---

## 8. Sky model

The external universe is represented by a **static equirectangular sky texture**.

- $\phi$ maps to the horizontal texture coordinate,
- $\theta$ maps to the vertical texture coordinate.

No frequency shift, Doppler effect, or radiative transfer is applied.  
The sampled sky color is used directly (up to gamma correction).

---

## 9. What is *not* implemented

The current implementation **does not include**:
- gravitational or Doppler redshift,
- intensity modulation via $k\cdot u$,
- emission from accretion disks or matter,
- adaptive error-controlled integration,
- explicit conservation checks (e.g. Carter constant).

These omissions are intentional and reflect the current scope of the project.

---

## 10. Scope

The physics model implemented here is intended to provide a **clean, coordinate-aware, strong-field ray-tracing framework** for Kerr spacetime, suitable for studying relativistic lensing and observer-dependent imaging effects.

It is a **first-principles geodesic simulation**, not an approximate visualization or post-Newtonian model.
