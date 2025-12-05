#include "kerr.hpp"
#include <cmath>

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

void rhs_null(const vec4& x, const vec4& p,
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
                  mat4& tet)
{
    // ---------------------------------------------------------
    // Metric at observer
    // ---------------------------------------------------------
    mat4 g = metric_cov(x_obs[0], x_obs[1], x_obs[2], x_obs[3], M, a);

    // ---------------------------------------------------------
    // 0) TIME LEG = normalized 4-velocity
    // ---------------------------------------------------------
    vec4 e0 = normalize_timelike(u_obs, g);

    // ---------------------------------------------------------
    // 1) SPATIAL SEEDS (coordinate basis)
    // ---------------------------------------------------------
    vec4 X = {0,1,0,0};   // radial (outward)
    vec4 Y = {0,0,1,0};   // polar
    vec4 Z = {0,0,0,1};   // azimuthal

    // ---------------------------------------------------------
    // Helper: project v orthogonal to e (metric)
    // ---------------------------------------------------------
    auto proj = [&](vec4 v, const vec4& e) {
        double c = ip(v, e, g);
        for(int mu=0; mu<4; ++mu)
            v[mu] -= c * e[mu];
        return v;
    };

    // ---------------------------------------------------------
    // 2) FORWARD VECTOR: inward radial direction
    // ---------------------------------------------------------
    vec4 e1 = X;        // start from outward radial
    e1 = proj(e1, e0);  // make it spatial in the rest frame
    e1 = normalize_spacelike(e1, g);

    // enforce INWARD orientation: want e1 opposite to X
    double s1 = ip(e1, X, g);
    if (s1 > 0.0) {
        for (int mu = 0; mu < 4; ++mu)
            e1[mu] = -e1[mu];
    }

    // ---------------------------------------------------------
    // 3) UP VECTOR: polar direction
    // ---------------------------------------------------------
    vec4 e2 = Y;
    e2 = proj(e2, e0);
    e2 = proj(e2, e1);
    e2 = normalize_spacelike(e2, g);

    // ---------------------------------------------------------
    // 4) RIGHT VECTOR: azimuthal direction
    // ---------------------------------------------------------
    vec4 e3 = Z;
    e3 = proj(e3, e0);
    e3 = proj(e3, e1);
    e3 = proj(e3, e2);
    e3 = normalize_spacelike(e3, g);

    // ---------------------------------------------------------
    // 5) APPLY CAMERA TILT (pure rotation in tetrad space)
    // ---------------------------------------------------------
    double eps = 0.0 * M_PI / 180.0;

    vec4 e1_new, e2_new;
    for(int mu=0; mu<4; ++mu) {
        e1_new[mu] = std::cos(eps)*e1[mu] + std::sin(eps)*e2[mu];
        e2_new[mu] = -std::sin(eps)*e1[mu] + std::cos(eps)*e2[mu];
    }

    e1 = normalize_spacelike(e1_new, g);
    e2 = proj(e2_new, e0);
    e2 = proj(e2, e1);
    e2 = normalize_spacelike(e2, g);

    // ---------------------------------------------------------
    // 6) STORE FINAL TETRAD
    // ---------------------------------------------------------
    for(int mu=0; mu<4; ++mu) {
        tet[0][mu] = e0[mu];  // time
        tet[1][mu] = e1[mu];  // forward (inward)
        tet[2][mu] = e2[mu];  // up
        tet[3][mu] = e3[mu];  // right
    }
}