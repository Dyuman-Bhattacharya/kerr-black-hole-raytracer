#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "kerr.hpp"
#include <string>
#include <omp.h>
#include <cmath>
#include <thread>

namespace py = pybind11;

// Compute horizon radius
static double horizon(double M, double a) {
    double discr = M*M - a*a;
    if (discr < 0) discr = 0;
    return M + std::sqrt(discr);
}

// Pixel → camera-frame direction (right, up, forward).
// z is simply the camera optical axis, NOT guaranteed to be radial inward/outward.


static void pixel_direction(int i, int j,
                            int width, int height,
                            double fov_deg,
                            double& nx, double& ny, double& nz)
{
    double aspect = double(width) / double(height);
    double fov = fov_deg * 3.1415926 / 180.0;

    // normalized device coordinates
    double x_ndc = (2.0 * (j + 0.5) / width  - 1.0) * aspect;
    double y_ndc = (1.0 - 2.0 * (i + 0.5) / height);
    double z     = 1.0 / std::tan(0.5 * fov);

    // forward = *INWARD* = +z in camera frame
    double len = std::sqrt(x_ndc*x_ndc + y_ndc*y_ndc + z*z);

    nx = x_ndc / len;     // right
    ny = y_ndc / len;     // up
    nz = z     / len;     // inward (correct)
}


static void build_p_from_pixel(int i, int j,
                               int width, int height,
                               double fov_deg,
                               const mat4& tetrad,
                               const mat4& g_cov,
                               vec4& p_cov_out)
{
    // k_tet = components in the camera tetrad frame (time, forward, up, right).
    // "forward" means +e1, which is physically the inward direction (because build_tetrad enforces that).
    // Do NOT interpret z as “outward”; it follows the tetrad, not camera geometry.


    double nx, ny, nz;
    pixel_direction(i, j, width, height, fov_deg, nx, ny, nz);

    // direction in CAMERA ORTHONORMAL FRAME
    // so nz attaches to +e1, not −e1
    vec4 k_tet = { -1.0,     // Set k_tet[0] = -1 only to break degeneracy; normalization is irrelevant for null rays.
                    nz,      
                    ny,      // up       (theta direction)
                    nx };    // right    (phi direction)

    // Transform to coordinate basis:  k^μ = e_(a)^μ k^(a)
    vec4 k_coord = {0,0,0,0};
    for (int a = 0; a < 4; ++a)
        for (int mu = 0; mu < 4; ++mu)
            k_coord[mu] += tetrad[a][mu] * k_tet[a];

    // Lower to covariant: p_μ = g_{μν} k^ν
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            s += g_cov[mu][nu] * k_coord[nu];
        p_cov_out[mu] = s;
    }
}

// ------------------------------------------------------------------
// Local helpers for timelike normalization in EF
// ------------------------------------------------------------------
static double ip_vec(const vec4& a, const vec4& b, const mat4& g)
{
    double s = 0.0;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            s += a[mu] * g[mu][nu] * b[nu];
    return s;
}

static void normalize_timelike_inplace(vec4& u, const mat4& g)
{
    double n = ip_vec(u, u, g);     // should be < 0 for timelike
    double s = std::sqrt(std::fabs(n));
    if (s == 0.0)
        throw std::runtime_error("normalize_timelike_inplace: null vector");
    for (int mu = 0; mu < 4; ++mu)
        u[mu] /= s;                 // g(u,u) ≈ -1 afterward (for n<0)
}

py::array_t<double> render_frame_full(
    py::array_t<double> x_obs_in, 
    py::array_t<double> u_obs_in,
    int width, int height,
    double fov_deg,
    double M, double a,
    double dl,
    int max_steps,
    double r_escape,
    py::array_t<double> sky_tex_in,
    double yaw_deg,
    double pitch_deg,
    double roll_deg)
{
    using std::sqrt;
    using std::pow;
    using std::max;

    omp_set_num_threads(std::thread::hardware_concurrency());

    const double pi        = 3.141592653589793;
    const double two_pi    = 2.0 * pi;
    const double inv_gamma = 1.0 / 2.2;

    // ------------------------------------------------------------
    // LOAD OBSERVER
    // ------------------------------------------------------------
    auto x0 = x_obs_in.unchecked<1>();
    auto u0 = u_obs_in.unchecked<1>();

    vec4 x_obs = {x0(0), x0(1), x0(2), x0(3)};
    vec4 u_obs = {u0(0), u0(1), u0(2), u0(3)};

    // Metric and tetrad at the camera position
    mat4 g_cov_cam = metric_cov(x_obs[0], x_obs[1], x_obs[2], x_obs[3], M, a);
    mat4 ginv_cam  = metric_contra(x_obs[0], x_obs[1], x_obs[2], x_obs[3], M, a);

    mat4 tetrad;
    build_tetrad(x_obs, u_obs, M, a, yaw_deg, pitch_deg, roll_deg, tetrad);
    //test_tetrad(x_obs, tetrad, M, a);

    // ------------------------------------------------------------
    // SKY TEXTURE
    // ------------------------------------------------------------
    auto sky_tex = sky_tex_in.unchecked<3>();
    py::ssize_t sky_h = sky_tex.shape(0);
    py::ssize_t sky_w = sky_tex.shape(1);

    // ------------------------------------------------------------
    // OUTPUT IMAGE
    // ------------------------------------------------------------
    py::array_t<double> rgb_final_arr({height, width, 3});
    auto rgb_final = rgb_final_arr.mutable_unchecked<3>();

    int    N      = width * height;
    double r_plus = horizon(M, a);
    double r_cap  = r_plus + 1e-2 * M;   // extremely tight gate
    bool camera_outside = (x_obs[1] > r_plus);

    // ------------------------------------------------------------
    // MAIN GEODESIC LOOP (one ray per pixel)
    // ------------------------------------------------------------
    #pragma omp parallel for
    for (int idx = 0; idx < N; ++idx)
    {
        int i = idx / width;
        int j = idx % width;

        // --------------------------------------------------------
        // Initial covariant momentum for this pixel
        // --------------------------------------------------------
        vec4 x = x_obs;
        double r_min_along = x[1];   // track minimum radius reached by this ray
        vec4 p;
        build_p_from_pixel(i, j, width, height, fov_deg,
                           tetrad, g_cov_cam, p);
        bool ever_outward = false;   // Did dr/dλ ever become positive?

        // --------------------------------------------------------
        // Integrate null geodesic
        // --------------------------------------------------------
        int    code      = 0;    // 0 = sky, 1 = horizon
        double theta_hit = x[2]; // will be overwritten on escape
        double phi_hit   = x[3];
        bool   finished  = false;

        for (int step = 0; step < max_steps; ++step)
        {
            double r_old = x[1];

            // Already in asymptotic region → treat current (θ,φ) as sky direction
            if (camera_outside && r_old > r_escape) {
                code      = 0;
                theta_hit = x[2];
                phi_hit   = x[3];
                finished  = true;
                break;
            }

            // Adaptive step size near hole / far out
            double dr_h   = r_old - r_plus;
            double dl_eff = dl;

            if (r_old > 80.0)      dl_eff *= 8.0;
            else if (r_old > 40.0) dl_eff *= 4.0;
            if (dr_h < 5.0) dl_eff *= 0.5;
            if (dr_h < 2.0) dl_eff *= 0.5;
            if (dr_h < 1.0) dl_eff *= 0.4;

            // RK4 step
            vec4 k1x, k2x, k3x, k4x;
            vec4 k1p, k2p, k3p, k4p;
            vec4 xm,  pm;

            geodesic_rhs(x, p, k1x, k1p, M, a);

            for (int mu = 0; mu < 4; ++mu) {
                xm[mu] = x[mu] + 0.5 * dl_eff * k1x[mu];
                pm[mu] = p[mu] + 0.5 * dl_eff * k1p[mu];
            }
            geodesic_rhs(xm, pm, k2x, k2p, M, a);

            for (int mu = 0; mu < 4; ++mu) {
                xm[mu] = x[mu] + 0.5 * dl_eff * k2x[mu];
                pm[mu] = p[mu] + 0.5 * dl_eff * k2p[mu];
            }
            geodesic_rhs(xm, pm, k3x, k3p, M, a);

            for (int mu = 0; mu < 4; ++mu) {
                xm[mu] = x[mu] + dl_eff * k3x[mu];
                pm[mu] = p[mu] + dl_eff * k3p[mu];
            }
            geodesic_rhs(xm, pm, k4x, k4p, M, a);

            for (int mu = 0; mu < 4; ++mu) {
                x[mu] += dl_eff * (k1x[mu] + 2.0*k2x[mu] + 2.0*k3x[mu] + k4x[mu]) / 6.0;
                p[mu] += dl_eff * (k1p[mu] + 2.0*k2p[mu] + 2.0*k3p[mu] + k4p[mu]) / 6.0;
            }

            // Compute dr/dλ using the inverse metric at the CURRENT ray position
            mat4 ginv_here = metric_contra(x[0], x[1], x[2], x[3], M, a);

            double dr_dl = 0.0;
            for (int nu = 0; nu < 4; ++nu)
                dr_dl += ginv_here[1][nu] * p[nu];

            if (dr_dl > 0.0)
                ever_outward = true;

            double r_new = x[1];
            if (r_new < r_min_along) r_min_along = r_new;

            // Horizon crossing: sign change in (r - r_plus)
            if ((r_old - r_plus) * (r_new - r_plus) <= 0.0 && r_old > r_plus) {
                code     = 1;
                finished = true;
                break;
            }

            // Just stepped beyond r_escape → record sky direction
            if (r_new > r_escape) {
                code      = 0;
                theta_hit = x[2];
                phi_hit   = x[3];
                finished  = true;
                break;
            }
        }

        if (camera_outside) {
            // Outside horizon: keep speckle fix
            if (r_min_along < r_cap)
                code = 1;
            else
                code = 0;
        } else {
            // Inside horizon: classify by turning point
            if (ever_outward)
                code = 0;   // sees external universe
            else
                code = 1;   // pure black
        }

        // --------------------------------------------------------
        // SAMPLE SKY TEXTURE (equirectangular) OR BLACK HOLE
        // --------------------------------------------------------
        double R = 0.0, G = 0.0, B = 0.0;

        if (code == 0 && sky_h > 0 && sky_w > 0)
        {
            // Normalize angles into [0,π] and [0,2π)
            double th = theta_hit;
            double ph = phi_hit;

            // Fold θ into [0,π]
            if (th < 0.0) th = -th;
            th = std::fmod(th, two_pi);
            if (th > pi) th = two_pi - th;

            // Wrap φ into [0,2π)
            ph = std::fmod(ph, two_pi);
            if (ph < 0.0) ph += two_pi;

            // Equirectangular mapping: φ -> x, θ -> y
            double u = ph / two_pi;  // [0,1)
            double v = th / pi;      // [0,1]

            double x_tex = u * (sky_w - 1);
            double y_tex = v * (sky_h - 1);

            int x0 = (int)std::floor(x_tex);
            int y0 = (int)std::floor(y_tex);

            double tx = x_tex - x0;
            double ty = y_tex - y0;

            if (x0 < 0) x0 = 0;
            if (x0 >= sky_w) x0 = (int)sky_w - 1;
            int x1 = x0 + 1;
            if (x1 >= sky_w) x1 = (int)sky_w - 1;

            if (y0 < 0) y0 = 0;
            if (y0 >= sky_h) y0 = (int)sky_h - 1;
            int y1 = y0 + 1;
            if (y1 >= sky_h) y1 = (int)sky_h - 1;

            // Bilinear sample
            for (int c = 0; c < 3; ++c)
            {
                double c00 = sky_tex(y0, x0, c);
                double c10 = sky_tex(y0, x1, c);
                double c01 = sky_tex(y1, x0, c);
                double c11 = sky_tex(y1, x1, c);

                double c0  = (1.0 - tx) * c00 + tx * c10;
                double c1  = (1.0 - tx) * c01 + tx * c11;
                double val = (1.0 - ty) * c0  + ty * c1;

                if (c == 0)      R = val;
                else if (c == 1) G = val;
                else             B = val;
            }
        }
        else
        {
            // If ray falls into hole or texture unavailable → return black pixel.
            R = G = B = 0.0;
        }

        // --------------------------------------------------------
        // GAMMA CORRECTION AND STORE
        // --------------------------------------------------------
        R = pow(max(R, 0.0), inv_gamma);
        G = pow(max(G, 0.0), inv_gamma);
        B = pow(max(B, 0.0), inv_gamma);

        rgb_final(i, j, 0) = R;
        rgb_final(i, j, 1) = G;
        rgb_final(i, j, 2) = B;
    }

    return rgb_final_arr;
}


// ------------------------------------------------------------------
// Integrate timelike geodesic in EF coordinates down to r_min
// Returns (tau_total, X[step_count,4], U[step_count,4])
// ------------------------------------------------------------------
py::tuple integrate_timelike_to_rmin(py::array_t<double> x0_in,
                                     py::array_t<double> u0_in,
                                     double M, double a,
                                     double d_tau,
                                     double r_min,
                                     int max_steps)
{
    // Load initial state
    auto x0v = x0_in.unchecked<1>();
    auto u0v = u0_in.unchecked<1>();

    vec4 x = { x0v(0), x0v(1), x0v(2), x0v(3) };   // (v,r,theta,phi)
    vec4 u = { u0v(0), u0v(1), u0v(2), u0v(3) };   // u^mu

    // Metric at start
    mat4 g_cov = metric_cov(x[0], x[1], x[2], x[3], M, a);
    mat4 g_inv = metric_contra(x[0], x[1], x[2], x[3], M, a);

    // Normalize u to timelike (g(u,u) = -1)
    normalize_timelike_inplace(u, g_cov);

    // Covariant momentum p_mu = g_{mu nu} u^nu
    vec4 p{};
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            s += g_cov[mu][nu] * u[nu];
        p[mu] = s;
    }

    // Storage for trajectory
    std::vector<vec4> X_traj;
    std::vector<vec4> U_traj;
    X_traj.reserve(max_steps + 1);
    U_traj.reserve(max_steps + 1);

    X_traj.push_back(x);
    U_traj.push_back(u);

    double tau = 0.0;

    for (int step = 0; step < max_steps; ++step)
    {
        double r = x[1];
        if (r <= r_min)
            break;

        // --- RK4 step for timelike geodesic using geodesic_rhs ---
        vec4 k1x, k2x, k3x, k4x;
        vec4 k1p, k2p, k3p, k4p;

        geodesic_rhs(x, p, k1x, k1p, M, a);

        vec4 xm, pm;

        for (int mu = 0; mu < 4; ++mu) {
            xm[mu] = x[mu] + 0.5 * d_tau * k1x[mu];
            pm[mu] = p[mu] + 0.5 * d_tau * k1p[mu];
        }
        geodesic_rhs(xm, pm, k2x, k2p, M, a);

        for (int mu = 0; mu < 4; ++mu) {
            xm[mu] = x[mu] + 0.5 * d_tau * k2x[mu];
            pm[mu] = p[mu] + 0.5 * d_tau * k2p[mu];
        }
        geodesic_rhs(xm, pm, k3x, k3p, M, a);

        for (int mu = 0; mu < 4; ++mu) {
            xm[mu] = x[mu] + d_tau * k3x[mu];
            pm[mu] = p[mu] + d_tau * k3p[mu];
        }
        geodesic_rhs(xm, pm, k4x, k4p, M, a);

        for (int mu = 0; mu < 4; ++mu) {
            x[mu] += d_tau * (k1x[mu] + 2*k2x[mu] + 2*k3x[mu] + k4x[mu]) / 6.0;
            p[mu] += d_tau * (k1p[mu] + 2*k2p[mu] + 2*k3p[mu] + k4p[mu]) / 6.0;
        }

        tau += d_tau;

        // --- Recompute metric at new point ---
        g_cov = metric_cov(x[0], x[1], x[2], x[3], M, a);
        g_inv = metric_contra(x[0], x[1], x[2], x[3], M, a);

        // --- Recover u^mu = g^{mu nu} p_nu ---
        for (int mu = 0; mu < 4; ++mu) {
            double s = 0.0;
            for (int nu = 0; nu < 4; ++nu)
                s += g_inv[mu][nu] * p[nu];
            u[mu] = s;
        }

        // --- Renormalize timelike & update p_mu ---
        normalize_timelike_inplace(u, g_cov);
        for (int mu = 0; mu < 4; ++mu) {
            double s = 0.0;
            for (int nu = 0; nu < 4; ++nu)
                s += g_cov[mu][nu] * u[nu];
            p[mu] = s;
        }

        X_traj.push_back(x);
        U_traj.push_back(u);
    }

    // Convert to numpy arrays
    py::ssize_t N = static_cast<py::ssize_t>(X_traj.size());

    py::array_t<double> X_arr({N, (py::ssize_t)4});
    py::array_t<double> U_arr({N, (py::ssize_t)4});

    auto X = X_arr.mutable_unchecked<2>();
    auto U = U_arr.mutable_unchecked<2>();

    for (py::ssize_t i = 0; i < N; ++i) {
        for (int mu = 0; mu < 4; ++mu) {
            X(i, mu) = X_traj[i][mu];
            U(i, mu) = U_traj[i][mu];
        }
    }

    return py::make_tuple(tau, X_arr, U_arr);
}


PYBIND11_MODULE(kerr_cpp, m) {
    m.doc() = "Fast Kerr geodesic integrator";

    m.def("metric_cov",
      [](double v, double r, double th, double ph, double M, double a) {
          mat4 g = metric_cov(v, r, th, ph, M, a);

          py::list outer;
          for (int i = 0; i < 4; ++i) {
              py::list inner;
              for (int j = 0; j < 4; ++j)
                  inner.append(g[i][j]);
              outer.append(inner);
          }

          return outer;
      },
      "Return covariant Kerr metric g_{mu nu} at (v,r,theta,phi)");

    m.def("render_frame_full", &render_frame_full,
        py::arg("x_obs"),
        py::arg("u_obs"),
        py::arg("width"),
        py::arg("height"),
        py::arg("fov_deg"),
        py::arg("M"),
        py::arg("a"),
        py::arg("dl"),
        py::arg("max_steps"),
        py::arg("r_escape"),
        py::arg("sky_texture"),
        py::arg("yaw_deg")   = 0.0,
        py::arg("pitch_deg") = 0.0,
        py::arg("roll_deg")  = 0.0);


    m.def("integrate_timelike_to_rmin", &integrate_timelike_to_rmin,
          py::arg("x0"),
          py::arg("u0"),
          py::arg("M"),
          py::arg("a"),
          py::arg("d_tau"),
          py::arg("r_min"),
          py::arg("max_steps"),
          "Integrate a timelike EF geodesic down to r_min; "
          "returns (tau_total, X[step,4], U[step,4]).");

}
