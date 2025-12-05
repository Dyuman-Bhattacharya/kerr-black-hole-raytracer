"""
Kerr Infall Ray Tracer
======================

A modular Python package for simulating null geodesics and visualizing
the optical appearance seen by an observer freely falling into a Kerr
black hole.

Modules:
    - metric:     Kerr metric and inverse in ingoing coordinates
    - geodesic:   Hamiltonian geodesic integrator (null & timelike)
    - observer:   Infalling observer worldline
    - tetrad:     Camera tetrad construction
    - scene:      Accretion disk geometry & emission
    - raytrace:   Per-ray integration logic
    - render_frame:  Frame rendering utilities
    - main:       Batch frame generation / movie creation
"""
