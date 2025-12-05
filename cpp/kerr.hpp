#pragma once
#include <array>

// Basic types
using vec4 = std::array<double, 4>;
using mat4 = std::array<std::array<double, 4>, 4>;

// Inverse Kerr metric g^{μν} in ingoing EF coordinates
mat4 metric_contra(double v, double r, double th, double ph,
                   double M, double a);

// Numerical partial derivatives ∂_α g^{μν}
std::array<mat4,4> metric_contra_derivs(double v, double r, double th, double ph,
                                        double M, double a);

// Hamiltonian RHS for null geodesics
void rhs_null(const vec4& x, const vec4& p,
              vec4& dx, vec4& dp,
              double M, double a);

// Covariant metric g_{μν}
mat4 metric_cov(double v, double r, double th, double ph, double M, double a);

// Build orthonormal tetrad e_(a)^μ at observer
// tetrad[a][mu] = e_(a)^μ, with e_(0)=u_obs,
// e_(1)=forward (inward), e_(2)=up, e_(3)=right.
void build_tetrad(const vec4& x_obs,
                  const vec4& u_obs,   // <-- NAME FIXED TO MATCH DEFINITION
                  double M, double a,
                  mat4& tetrad);