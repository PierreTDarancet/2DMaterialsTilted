# Tilted Supercell Generator for 2D Materials

A Python tool for generating low-symmetry, orthorhombic supercells of 2D materials (specifically MoS₂) tilted under periodic boundary conditions.

## 📌 Problem Overview

In Density Functional Theory (DFT) calculations involving finite electric fields (e.g., Berry phase, sawtooth potentials) or specific dislocation geometries, it is often necessary to break the standard alignment of the 2D material with the simulation cell axes while maintaining **Periodic Boundary Conditions (PBC)** and an **orthorhombic** cell shape.

This tool solves the geometry problem of generating a supercell where:
1.  The simulation box remains strictly orthorhombic ($\alpha=\beta=\gamma=90^\circ$).
2.  The 2D material sheet is tilted by an angle $\theta$ relative to the cell axes.
3.  The perpendicular distance (vacuum spacing) between periodic images is exactly $h$.
4.  Internal bond lengths and angles (e.g., the Mo-S trigonal prisms) are preserved via rigid-body rotation.

---

## 📐 Geometric Derivation

### 1. The "Hypotenuse" Constraint
Consider a 2D material with lattice constant $a_0$. We wish to create a supercell where the material runs along the diagonal of a new box. For the PBCs to connect perfectly, the length of the diagonal must match an integer supercell size $N$ of the material:

$$L_{\text{sheet}} = \sqrt{L_{\text{planar}}^2 + L_z^2} = N \cdot a_0$$

### 2. The Spacing Constraint
The perpendicular distance $h$ between the periodic images (the sheet and its copy) is determined by the area of the cell divided by the length of the diagonal:

$$h = \frac{L_{\text{planar}} L_z}{N \cdot a_0}$$

### 3. The Solution
Solving this system for the box dimensions ($L_{\text{planar}}, L_z$) yields a quadratic equation. A real physical solution exists only if the **Stability Condition** is met:

$$N \cdot a_0 \geq 2h$$

If this condition is met, the box dimensions are:

$$L_{\text{planar}} = \sqrt{\frac{(N a_0)^2 + \sqrt{(N a_0)^4 - 4(N a_0 h)^2}}{2}}$$

$$L_z = \frac{N a_0 h}{L_{\text{planar}}}$$

---

## ⚙️ Methodology

### Phase Initialization (MoS₂ 2H)
To prevent distortions, the code initializes the **Centered Rectangular Approximant** of the hexagonal cell:
* **Lattice:** $a_{\text{rect}} = a_{\text{hex}}$, $b_{\text{rect}} = a_{\text{hex}}\sqrt{3}$
* **Basis:**
    * Mo at $(0,0)$ and $(0.5, 0.5)$
    * S at $(0, 1/3)$ and $(0.5, 5/6)$ relative to Mo.
    * *This ensures correct trigonal prismatic coordination.*

### Rigid Body Rotation
To preserve finite thickness (e.g., the S-Mo-S sandwich height), we perform a rigid Cartesian rotation rather than coordinate scaling.

For a rotation around the Y-axis (tilting the sheet in the XZ plane), we apply the rotation matrix $\mathbf{R}_y(\theta)$ to every atom's Cartesian vector $(x,y,z)$ relative to the slab center:

$$
\mathbf{r}' = 
\begin{pmatrix}
\cos\theta & 0 & -\sin\theta \\
0 & 1 & 0 \\
\sin\theta & 0 & \cos\theta
\end{pmatrix}
\mathbf{r}
$$

---

## 🚀 Usage

### Dependencies
* Python 3.x
* `numpy`
* `pymatgen`

### Running the Code
The tool automatically calculates the minimum supercell size $N$ required to satisfy your requested vacuum spacing $h$.

```python
# Example: Generate a cell with 15.0 Angstrom spacing
# rotating around the b-axis (tilting the 'a' vector).

generate_tilted_cell(target_h=15.0, rotation_axis='b')
