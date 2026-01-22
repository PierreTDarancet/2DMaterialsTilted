# @title Rigid-Body Tilted Supercell Generator (Corrected Geometry)
# @markdown Run this to generate a strictly orthorhombic tilted cell with preserved bond angles.

import numpy as np
import os
import subprocess
import sys

# --- 1. SETUP ---
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

try:
    import pymatgen.core
except ImportError:
    print("Installing pymatgen...")
    install("pymatgen")

from pymatgen.core import Structure, Lattice
from google.colab import files

# --- 2. GENERATION LOGIC ---

def get_rectangular_mos2(a_lattice=3.18, thickness=3.12, vacuum=20.0):
    """
    Generates the standard 6-atom Rectangular MoS2 cell (Orthorhombic).
    a_rect = a
    b_rect = a * sqrt(3)
    """
    b_lattice = a_lattice * np.sqrt(3)
    
    # Orthorhombic Lattice
    lattice = Lattice.from_parameters(a_lattice, b_lattice, vacuum, 90, 90, 90)
    
    # Fractional coordinates for 2H-MoS2 in rectangular setting
    # Mo at (0, 0) and (0.5, 0.5) roughly
    # This is the standard 'centered rectangular' representation
    
    # Thickness in fractional z
    z_delta = thickness / vacuum / 2.0
    z_mo = 0.5
    z_s_top = 0.5 + z_delta
    z_s_bot = 0.5 - z_delta
    
    species = ["Mo", "Mo", "S", "S", "S", "S"]
    coords = [
        [0.0, 0.0, z_mo],          # Mo1
        [0.5, 0.5, z_mo],          # Mo2
        [0.0, 0.0, z_s_bot],       # S1 Bottom
        [0.0, 0.0, z_s_top],       # S1 Top
        [0.5, 0.5, z_s_bot],       # S2 Bottom
        [0.5, 0.5, z_s_top],       # S2 Top
    ]
    
    return Structure(lattice, species, coords)

def solve_box_geometry(length_sheet, h_spacing):
    """Calculates new box dimensions and tilt angle theta."""
    if length_sheet < 2 * h_spacing:
        raise ValueError(f"Geometry impossible: Sheet length ({length_sheet:.2f}) < 2*h ({2*h_spacing}). Increase N.")
        
    # Quadratic solution for box dimensions
    term = np.sqrt(length_sheet**4 - 4 * (length_sheet * h_spacing)**2)
    L_planar = np.sqrt((length_sheet**2 + term) / 2)
    L_z = (length_sheet * h_spacing) / L_planar
    theta = np.arctan(L_z / L_planar) # Angle in radians
    
    return L_planar, L_z, theta

def create_tilted_cell(N, h_spacing, a_param=3.18):
    # 1. Start with the perfect Rectangular Cell (6 atoms)
    struct = get_rectangular_mos2(a_lattice=a_param)
    
    # 2. Create the "Strip" (Supercell along X)
    # We always rotate around Y for this recipe (tilting X into Z).
    # This is mathematically equivalent to rotating around X if we just swapped labels,
    # so we standardize on Y-rotation to keep the math clean.
    struct.make_supercell([N, 1, 1])
    
    # Get Cartesian Coords (Orthogonal basis guarantees these are simple x,y,z)
    cart_coords = struct.cart_coords
    
    # Center the strip at z=0 relative to the slab center
    # This ensures rotation happens "around the middle" of the slab
    z_mean = np.mean(cart_coords[:, 2])
    cart_coords[:, 2] -= z_mean
    
    # 3. Solve Geometry
    total_length = struct.lattice.a # This is N * 3.18
    L_planar, L_z, theta = solve_box_geometry(total_length, h_spacing)
    
    # 4. Rigid Body Rotation (Cartesian)
    # Rotation matrix for rotating around Y-axis (X -> Z)
    c, s = np.cos(theta), np.sin(theta)
    
    # R_y = [[c,  0, -s],
    #        [0,  1,  0],
    #        [s,  0,  c]]
    
    rotation_matrix = np.array([
        [c, 0, -s],
        [0, 1, 0],
        [s, 0, c]
    ])
    
    # Apply rotation to all vectors
    # Transpose needed because numpy dots last axis
    rotated_coords = np.dot(cart_coords, rotation_matrix.T)
    
    # 5. Define New Lattice
    # The new lattice is Orthorhombic
    old_b = struct.lattice.b
    new_lattice = Lattice.from_parameters(L_planar, old_b, L_z, 90, 90, 90)
    
    # 6. Wrap Coordinates
    # Convert to fractional coordinates in the NEW box to handle wrapping easily
    new_struct_temp = Structure(new_lattice, struct.species, rotated_coords, coords_are_cartesian=True)
    
    # This step automatically wraps atoms back into the [0,1] box
    # effectively handling the periodic boundary conditions for us
    final_struct = new_struct_temp.copy()
    
    return final_struct, theta

# --- 3. EXECUTION ---

# Parameters
N_SUPERCELL = 10     # Number of rectangular cells (Length will be N * 3.18)
H_SPACING = 15.0     # Vacuum spacing
A_PARAM = 3.18       # Lattice constant of MoS2

try:
    # Generate
    tilted_struct, theta_rad = create_tilted_cell(N_SUPERCELL, H_SPACING, A_PARAM)
    
    # Write output
    output_file = "POSCAR_tilted_ortho.vasp"
    tilted_struct.to(fmt="poscar", filename=output_file)
    
    print(f"SUCCESS: {output_file} generated.")
    print(f"Structure: {tilted_struct.formula}")
    print(f"Tilt Angle: {np.degrees(theta_rad):.2f} degrees")
    print(f"New Lattice: {tilted_struct.lattice.abc}")
    
    # Download
    files.download(output_file)

except Exception as e:
    print(f"ERROR: {e}")
