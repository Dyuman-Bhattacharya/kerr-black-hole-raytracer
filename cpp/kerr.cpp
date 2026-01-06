#include "kerr.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <pybind11/pybind11.h>
namespace py = pybind11;


mat4 metric_contra(double v, double r, double th, double ph, double M, double a) {
    (void)v; (void)ph;  // not used explicitly

    mat4 ginv{};  // zero-initialize

    double sinth = std::sin(th);
    double sin2  = sinth * sinth;
    double costh = std::cos(th);
    double cos2  = costh * costh;

    double Sigma = r*r + a*a * cos2;
    double Delta = r*r - 2.0*M*r + a*a;

    // Correct ingoing EF inverse metric g^{μν} in (v, r, θ, φ)
    // g^{vv}
    ginv[0][0] = (a*a * sin2) / Sigma;

    // g^{vr} = g^{rv}
    ginv[0][1] = (a*a + r*r) / Sigma;
    ginv[1][0] = ginv[0][1];

    // g^{rr}
    ginv[1][1] = Delta / Sigma;

    // g^{vφ} = g^{φv}
    ginv[0][3] = a / Sigma;
    ginv[3][0] = ginv[0][3];

    // g^{rφ} = g^{φr}
    ginv[1][3] = a / Sigma;
    ginv[3][1] = ginv[1][3];

    // g^{θθ}
    ginv[2][2] = 1.0 / Sigma;

    // g^{φφ}
    ginv[3][3] = 1.0 / (Sigma * sin2);

    return ginv;
}

std::array<mat4,4> metric_contra_derivs(
    double v, double r, double th, double ph, double M, double a)
{
    std::array<mat4,4> dginv{};

    double sinth = std::sin(th);
    double costh = std::cos(th);
    double sin2  = sinth*sinth;
    double cos2  = costh*costh;

    double Sigma = r*r + a*a*cos2;
    double Delta = r*r - 2*M*r + a*a;

    double dSigma_dr  = 2*r;
    double dSigma_dth = -2*a*a*sinth*costh;

    double dDelta_dr  = 2*r - 2*M;

    double Sigma2 = Sigma*Sigma;

    // ----------------------------
    // ∂_r g^{μν}
    // ----------------------------
    mat4& dr = dginv[1];

    dr[0][0] = -(a*a*sin2)*dSigma_dr / Sigma2;

    dr[0][1] = dr[1][0] =
        ( (2*r)*Sigma - (r*r + a*a)*dSigma_dr ) / Sigma2;

    dr[1][1] = ( dDelta_dr*Sigma - Delta*dSigma_dr ) / Sigma2;

    dr[0][3] = dr[3][0] = -(a*dSigma_dr) / Sigma2;

    dr[1][3] = dr[3][1] = -(a*dSigma_dr) / Sigma2;

    dr[2][2] = -(dSigma_dr) / Sigma2;

    dr[3][3] = -(dSigma_dr) / (Sigma2 * sin2);
    

    // ----------------------------
    // ∂_θ g^{μν}
    // ----------------------------
    mat4& dt = dginv[2];

    dt[0][0] =
        -(a*a*sin2)*dSigma_dth / Sigma2
        + (a*a * 2*sinth*costh) / Sigma;

    dt[0][1] = dt[1][0] =
        -(r*r + a*a)*dSigma_dth / Sigma2;

    dt[1][1] =
        -(Delta * dSigma_dth) / Sigma2;

    dt[0][3] = dt[3][0] =
        -(a*dSigma_dth) / Sigma2;

    dt[1][3] = dt[3][1] =
        -(a*dSigma_dth) / Sigma2;

    dt[2][2] =
        -(dSigma_dth) / Sigma2;
        
    dt[3][3] =
        -(dSigma_dth) / (Sigma2 * sin2)
        - 2.0 * costh / (Sigma * sinth * sin2);

    return dginv;
}

void geodesic_rhs(const vec4& x, const vec4& p,
              vec4& dx, vec4& dp,
              double M, double a) {
    double v  = x[0];
    double r  = x[1];
    double th = x[2];
    double ph = x[3];

    mat4 ginv = metric_contra(v, r, th, ph, M, a);
    auto dg   = metric_contra_derivs(v, r, th, ph, M, a);

    // dx^μ = g^{μν} p_ν
    for (int mu = 0; mu < 4; ++mu) {
        dx[mu] = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            dx[mu] += ginv[mu][nu] * p[nu];
    }

    // dp_μ = -1/2 ∂_μ g^{αβ} p_α p_β
    for (int mu = 0; mu < 4; ++mu) {
        double acc = 0.0;
        for (int a_i = 0; a_i < 4; ++a_i)
            for (int b_i = 0; b_i < 4; ++b_i)
                acc += dg[mu][a_i][b_i] * p[a_i] * p[b_i];

        dp[mu] = -0.5 * acc;
    }
}

mat4 metric_cov(double v, double r, double th, double ph, double M, double a)
{
    (void)v; (void)ph; // not used explicitly

    mat4 g{};

    double sinth = std::sin(th);
    double sin2  = sinth * sinth;
    double costh = std::cos(th);
    double cos2  = costh * costh;

    double Sigma = r*r + a*a * cos2;
    double f     = 2.0 * M * r / Sigma;

    // Indices: 0=v, 1=r, 2=theta, 3=phi~
    // g_vv
    g[0][0] = -(1.0 - f);

    // g_vr = g_rv
    g[0][1] = 1.0;
    g[1][0] = 1.0;

    // g_vphi = g_phiv
    g[0][3] = -2.0 * a * M * r * sin2 / Sigma;
    g[3][0] = g[0][3];

    // g_rr = 0 (already zero)

    // g_rphi = g_phir
    g[1][3] = -a * sin2;
    g[3][1] = g[1][3];

    // g_thetatheta
    g[2][2] = Sigma;

    // g_phiphi
    g[3][3] = (r*r + a*a + 2.0 * M * a*a * r * sin2 / Sigma) * sin2;

    return g;
}

static double ip(const vec4& a, const vec4& b, const mat4& g)
{
    double s = 0.0;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            s += a[mu] * g[mu][nu] * b[nu];
    return s;
}

static vec4 normalize_timelike(const vec4& u, const mat4& g)
{
    double n = ip(u, u, g);  // should be negative
    double s = std::sqrt(std::fabs(n));
    vec4 out{};
    for (int mu = 0; mu < 4; ++mu) out[mu] = u[mu] / s;
    return out;
}

static vec4 normalize_spacelike(const vec4& v, const mat4& g)
{
    double n = ip(v, v, g);  // should be positive
    double s = std::sqrt(std::fabs(n));
    vec4 out{};
    for (int mu = 0; mu < 4; ++mu) out[mu] = v[mu] / s;
    return out;
}

// project v into rest space of u: v -> v + (u·v) u   (since u·u = -1)
static vec4 project_orthogonal(const vec4& v, const vec4& u, const mat4& g)
{
    double uv = ip(u, v, g);
    vec4 out{};
    for (int mu = 0; mu < 4; ++mu)
        out[mu] = v[mu] + uv * u[mu];
    return out;
}

float M_PI=3.1415926;

//==============================================================
//  TETRAD BUILDER (EF–Kerr) — CLEAN, CORRECT, ORTHONORMAL
//==============================================================
void build_tetrad(const vec4& x_obs,
                  const vec4& u_obs,
                  double M, double a,
                  double yaw_deg, double pitch_deg, double roll_deg,
                  mat4& tet)
{
    // Metric and inverse at observer position
    mat4 g     = metric_cov  (x_obs[0], x_obs[1], x_obs[2], x_obs[3], M, a);
    mat4 g_inv = metric_contra(x_obs[0], x_obs[1], x_obs[2], x_obs[3], M, a);

    // Time leg: normalized 4-velocity (future-directed timelike)
    vec4 e0 = normalize_timelike(u_obs, g);  // g(e0,e0) = -1

    // Lorentzian Gram–Schmidt projector:
    // v -> v - [(v·e)/(e·e)] e   (works for timelike or spacelike e)
    auto proj = [&](vec4 v, const vec4& e) -> vec4 {
        double c = ip(v, e, g);
        double n = ip(e, e, g);
        if (std::fabs(n) < 1e-14)
            return v;  // don't try to project on (numerically) null vector
        for (int mu = 0; mu < 4; ++mu)
            v[mu] -= (c / n) * e[mu];
        return v;
    };

    // ---------------------------------------------------------------------
    // 1. Forward leg e1: roughly radial direction, orthogonal to e0
    // ---------------------------------------------------------------------
    vec4 e1_seed = {0.0, 1.0, 0.0, 0.0};  // coordinate r-direction
    vec4 e1 = proj(e1_seed, e0);          // remove time component
    e1 = normalize_spacelike(e1, g);      // g(e1,e1) = +1

    // ---------------------------------------------------------------------
    // 2. Up leg e2: θ-direction, orthogonal to e0 and e1
    // ---------------------------------------------------------------------
    vec4 e2_seed = {0.0, 0.0, 1.0, 0.0};
    vec4 e2 = proj(e2_seed, e0);
    e2 = proj(e2, e1);
    e2 = normalize_spacelike(e2, g);

    // ---------------------------------------------------------------------
    // 3. Right leg e3: φ-direction, orthogonal to e0, e1, e2
    // ---------------------------------------------------------------------
    vec4 e3_seed = {0.0, 0.0, 0.0, 1.0};
    vec4 e3 = proj(e3_seed, e0);
    e3 = proj(e3, e1);
    e3 = proj(e3, e2);
    e3 = normalize_spacelike(e3, g);

    // ---------------------------------------------------------------------
    // 4. Ensure "forward" (e1) is actually inward: dr/dλ < 0 for -e0+e1
    // ---------------------------------------------------------------------
    vec4 k_coord{};
    for (int mu = 0; mu < 4; ++mu)
        k_coord[mu] = -e0[mu] + e1[mu];   // null direction candidate

    // p_mu = g_{μν} k^ν
    vec4 p_test{};
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            s += g[mu][nu] * k_coord[nu];
        p_test[mu] = s;
    }

    // dr/dλ = k^r = g^{rν} p_ν   (index 1 is r)
    double dr = 0.0;
    for (int nu = 0; nu < 4; ++nu)
        dr += g_inv[1][nu] * p_test[nu];

    // If this null direction goes outward (dr > 0), flip spatial triad
    if (dr > 0.0) {
        for (int mu = 0; mu < 4; ++mu) {
            e1[mu] = -e1[mu];
            e2[mu] = -e2[mu];
            e3[mu] = -e3[mu];
        }
    }

    // ---------------------------------------------------------------------
    // 5. Apply yaw / pitch / roll to spatial triad (e1,e2,e3)
    //    Uses existing global M_PI from your file.
    // ---------------------------------------------------------------------
    double cy = std::cos(yaw_deg   * M_PI / 180.0);
    double sy = std::sin(yaw_deg   * M_PI / 180.0);
    double cp = std::cos(pitch_deg * M_PI / 180.0);
    double sp = std::sin(pitch_deg * M_PI / 180.0);
    double cr = std::cos(roll_deg  * M_PI / 180.0);
    double sr = std::sin(roll_deg  * M_PI / 180.0);

    // Same rotation convention you already had
    double R[3][3] = {
        { cy*cr + sy*sp*sr,   sr*cp,   -sy*cr + cy*sp*sr },
        { -cy*sr + sy*sp*cr,  cr*cp,    sr*sy + cy*sp*cr },
        { sy*cp,              -sp,      cy*cp            }
    };

    vec4 e1_r{}, e2_r{}, e3_r{};
    for (int mu = 0; mu < 4; ++mu) {
        double v1 = e1[mu];
        double v2 = e2[mu];
        double v3 = e3[mu];

        e1_r[mu] = R[0][0]*v1 + R[0][1]*v2 + R[0][2]*v3;
        e2_r[mu] = R[1][0]*v1 + R[1][1]*v2 + R[1][2]*v3;
        e3_r[mu] = R[2][0]*v1 + R[2][1]*v2 + R[2][2]*v3;
    }

    // Re-orthonormalize after rotation, just to kill numerical drift
    e1 = normalize_spacelike(e1_r, g);
    e2 = proj(e2_r, e1);
    e2 = normalize_spacelike(e2, g);
    e3 = proj(e3_r, e1);
    e3 = proj(e3, e2);
    e3 = normalize_spacelike(e3, g);

    // ---------------------------------------------------------------------
    // 6. Store tetrad as rows: tetrad[a][mu] = e_(a)^μ
    // ---------------------------------------------------------------------
    for (int mu = 0; mu < 4; ++mu) {
        tet[0][mu] = e0[mu];  // time
        tet[1][mu] = e1[mu];  // forward (inward)
        tet[2][mu] = e2[mu];  // up
        tet[3][mu] = e3[mu];  // right
    }
}

//==============================================================
//   TETRAD DIAGNOSTIC — RUN THIS AFTER build_tetrad()
//==============================================================

void test_tetrad(const vec4& x,
                 const mat4& tetrad,
                 double M, double a)
{
    mat4 g     = metric_cov  (x[0], x[1], x[2], x[3], M, a);
    mat4 g_inv = metric_contra(x[0], x[1], x[2], x[3], M, a);

    // Compute G_ab = g_{μν} e_a^μ e_b^ν
    double G[4][4];
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            double s = 0.0;
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    s += tetrad[a][mu] * g[mu][nu] * tetrad[b][nu];
            G[a][b] = s;
        }
    }

    bool ok = true;

    // --- Check time leg ---
    if (std::fabs(G[0][0] + 1.0) > 1e-6) {
        py::print("[TETRAD ERROR] g(e0,e0) =", G[0][0], "(expected -1)");
        ok = false;
    }

    // --- Check spatial legs ---
    for (int i = 1; i < 4; ++i) {
        if (std::fabs(G[i][i] - 1.0) > 1e-6) {
            py::print("[TETRAD ERROR] g(e", i, ",e", i, ") =", G[i][i], "(expected +1)");
            ok = false;
        }
    }

    // --- Check cross terms ---
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            if (a == b) continue;
            if (std::fabs(G[a][b]) > 1e-6) {
                py::print("[TETRAD ERROR] g(e", a, ",e", b, ") =", G[a][b], "(expected 0)");
                ok = false;
            }
        }
    }

    // --- Check inward-pointing direction ---
    vec4 e0{}, e1{};
    for (int mu = 0; mu < 4; ++mu) {
        e0[mu] = tetrad[0][mu];
        e1[mu] = tetrad[1][mu];
    }

    vec4 k{};  // null direction candidate
    for (int mu = 0; mu < 4; ++mu)
        k[mu] = -e0[mu] + e1[mu];

    vec4 p{};
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            s += g[mu][nu] * k[nu];
        p[mu] = s;
    }

    double dr = 0.0;
    for (int nu = 0; nu < 4; ++nu)
        dr += g_inv[1][nu] * p[nu];

    if (dr >= 0.0) {
        py::print("[TETRAD ERROR] forward leg is NOT inward: dr/dλ =", dr, "(expected < 0)");
        ok = false;
    }

    // --- Final report ---
    if (ok) {
        py::print("[TETRAD OK] Orthonormal and inward.");
    }
}
