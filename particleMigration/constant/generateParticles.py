#!/usr/bin/env python3
"""
Generate an OpenFOAM kinematicCloudPositions file.

Particles are placed one per Cartesian sub-cell, with a random perturbation
that:
    1. Keeps each particle inside its own sub-cell.
    2. Prevents overlap with particles in neighboring cells.

The script also computes the particle volume fraction.

Author: ChatGPT
"""

import numpy as np

# =============================================================================
# User parameters
# =============================================================================

# Box limits (m)
xmin, xmax = 0.0, 785e-6
ymin, ymax = 0.0, 785e-6
zmin, zmax = 0.0, 785e-6

# Cell spacing (m)
deltax = xmax/10
deltay = xmax/10
deltaz = xmax/10

# Particle diameter (m)
k = 2*0.0365 #0.0365
diameter = xmax*k

# Random seed (None for different realizations every run)
seed = 12345

# Output file
outfile = "kinematicCloudPositions"
outfiled = "d.H"

# =============================================================================
# Initialize random number generator
# =============================================================================

rng = np.random.default_rng(seed)

# =============================================================================
# Derived quantities
# =============================================================================

Lx = xmax - xmin
Ly = ymax - ymin
Lz = zmax - zmin

nx = int(np.floor(Lx / deltax))
ny = int(np.floor(Ly / deltay))
nz = int(np.floor(Lz / deltaz))

print("Grid:")
print(f"  nx = {nx}")
print(f"  ny = {ny}")
print(f"  nz = {nz}")

# Check feasibility
if diameter > min(deltax, deltay, deltaz):
    raise ValueError(
        "Particle diameter is larger than the cell spacing.\n"
        "Particles in neighboring cells would necessarily overlap."
    )

# =============================================================================
# Maximum admissible perturbation
# =============================================================================

# Limited by cell boundaries
rx_cell = deltax / 2
ry_cell = deltay / 2
rz_cell = deltaz / 2

# Limited by overlap constraint
rx_overlap = (deltax - diameter) / 2
ry_overlap = (deltay - diameter) / 2
rz_overlap = (deltaz - diameter) / 2

# Final allowable perturbation
rx = min(rx_cell, rx_overlap)
ry = min(ry_cell, ry_overlap)
rz = min(rz_cell, rz_overlap)

print("\nMaximum random perturbation:")
print(f"  x : ±{rx:.6e} m")
print(f"  y : ±{ry:.6e} m")
print(f"  z : ±{rz:.6e} m")

# =============================================================================
# Generate particle positions
# =============================================================================

positions = []

for k in range(nz):

    zc = zmin + (k + 0.5) * deltaz

    for j in range(ny):

        yc = ymin + (j + 0.5) * deltay

        for i in range(nx):

            xc = xmin + (i + 0.5) * deltax

            x = xc + rng.uniform(-rx, rx)
            y = yc + rng.uniform(-ry, ry)
            z = zc + rng.uniform(-rz, rz)

            positions.append((x, y, z))

nParticles = len(positions)

# =============================================================================
# Volume fraction
# =============================================================================

particle_volume = np.pi * diameter**3 / 6.0
total_particle_volume = nParticles * particle_volume
box_volume = Lx * Ly * Lz

volume_fraction = total_particle_volume / box_volume

inter_distance = 2*diameter

# =============================================================================
# Write OpenFOAM file
# =============================================================================

with open(outfile, "w") as f:

    f.write(
"""FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant";
    object      kinematicCloudPositions;
}

"""
    )

    f.write(f"{nParticles}\n")
    f.write("(\n")

    for p in positions:
        f.write(f"({p[0]:.12g} {p[1]:.12g} {p[2]:.12g})\n")

    f.write(")\n")

with open(outfiled, "w") as f:
    f.write(f"dp {diameter};\n")
    f.write(f"idist {inter_distance};\n")

# =============================================================================
# Summary
# =============================================================================

print("\n==========================================")
print("OpenFOAM particle cloud generated")
print("==========================================")
print(f"Box dimensions          : {Lx:.6e} x {Ly:.6e} x {Lz:.6e} m")
print(f"Cell spacing            : {deltax:.6e}, {deltay:.6e}, {deltaz:.6e} m")
print(f"Particle diameter       : {diameter:.6e} m")
print(f"Number of particles     : {nParticles}")
print(f"Particle volume         : {particle_volume:.6e} m³")
print(f"Total particle volume   : {total_particle_volume:.6e} m³")
print(f"Box volume              : {box_volume:.6e} m³")
print(f"Volume fraction         : {volume_fraction:.6f}")
print("==========================================")
