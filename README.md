# Vicsek-like ABM of cell aggregation in cadherin-linked synthetic switch

Link to manuscript: To be added

## Title
Synthetic symmetry breaking and programmable multicellular structure formation
## Authors
*Noreen Wauford, *Akshay Patel, *Jesse Tordoff, *Casper Enghuus, Andrew Jin<sup>++</sup>, Jack Toppen, Melissa L. Kemp, Ron Weiss

<sup>++</sup>: Code author

# Model Description
Cells are modeled as rigid spheres of three cell types: blue, red, and yellow cells all with equal cell diameter σ. The three cell types vary in adhesion strength, and are in ascending order: blue, red, and yellow, with yellow-yellow adhesion the highest. All cells are initialized as blue cells with a fixed probability of transitioning to a cadherin expressing (red, yellow) state determined by the final fluorescence intensity expression of the corresponding Dox/ABA dosage experiment, p_total cells committed. The fixed transition probability at each time step (x) is back calculated from the final fluorescence intensity expression:

$$(1-x)^{72h}=(1-p_{total cells committed})$$

When a cell is indicated to transition to a cadherin expressing state, the probability of choosing the red or yellow states is also determined by the fluorescence intensity expression. The probability of choosing the one of the states is equal to the experimental intensity expression ratio of the respective cells over the total non-blue intensity expression. 

Cell locations are stored in a NumPy array and used for local cell interaction calculations. Cell locations are updated each time step (10 simulated seconds), at a fixed velocity of 0.3 * σ (cell diameter). The direction of each cell’s velocity vector is determined by the normalized sum of all forces acting upon each cell. These forces are Brownian motion, a gravity field, and the local cell-cell interactions of homotypic and heterotypic adhesion and steric repulsion. 

All cells experience Brownian motion and a gravity field designed to mimic cell behavior in an ultra-low binding U-shaped bottom well. To represent the curvature of the well, the magnitude of the gravity field is dependent on the xy planar distance $(|r_{i,xy}|)$ from the center of the simulation space. $U_g$ represents the strength of the force and R as the radius of the well.

$$f_i=-U_g r_i  (|r_{i,xy}|)/R$$

The presence of the pairwise forces is dependent on the distance $r_{ij}$ between two cells. At $|r_{ij}<1|$, two cells are overlapping, and each cell experiences a strong $(U_{hc}=10^4)$ steric repulsion force to remove overlap of the rigid bodies.

$$f_{ij}=-U_{hc}  r_{ij}/|r_{ij}|$$ 

At $|1< r_{ij}<1.6*σ|$, cells experience attraction due to cell-cell adhesion. This is described by

$$f_{ij}=U_{ij} (r_{ij}-r_e )  r_{ij}/|r_{ij}|$$

Where $U_{ij}$ is the pairwise cell type specific adhesion force. There are 6 cell-type specific adhesion parameters, of which 3 are for the homotypic $(U_{bb},U_{rr},U_{yy})$ adhesion forces, and 3 are for the heterotypic $(U_{br},U_{by},U_{ry})$ adhesion forces. The equilibrium distance of intercellular interaction is denoted as $r_e$. The isotropic noise of the attraction force is represented with a Gaussian noise vector $ξ_{ij}$ with magnitude α. At  $|r_{ij}>1.6*σ|$, the cells are sufficiently far and do not interact. As cells sufficiently distant from each other do not influence each other, local interactions are identified through a fixed radius search of 1.6 * σ for nearby cell bodies and then subsequently divided into adhesion and repulsion forces. 

# Implementation
We implemented the ABM using the PythonABM library developed by author Jack Toppen. PythonABM is an installable library that provides a framework for building efficient ABMs in Python. 

PythonABM is hosted on Python’s package index, PyPI, such that the library can be installed through the command-line. For more information, PythonABM’s documentation describes the complete installation process and provides an example script for starting to use the framework. PythonABM is available at [here](https://github.com/kemplab/pythonabm)  or on [PyPI](https://pypi.org/project/pythonabm/). 

# How to run
To run a single simulation, the easiest way is to run 
```
python model.py rr yy ry dox aba tf subts
```
RR, YY, RY specifies the respective cell-cell adhesion parameters. <br>
dox, aba are the final ratios of red and yellow cells, respectively <br>
tf is the final number of frames <br>
subts is the number of time steps per frame. <br>

This command will create a simulation in the current directory. The output files will be:
1. A folder containing csv files of the locations of the cells at each frame
2. A folder contains png files of the cells at each frame
3. A csv file containing the duration of each timestep
4. A mov file containing the resultant simulation
5. A .pckl file which allows you to resume the simulation if interrupted.

To recreate simulations as seen in Fig 4 of the article, set RR=30, YY=40, RY=1. The dox and aba ratios vary per square in the matrix. To run all the simulations found in the replicates found in Fig 4, the following method can be run instead.
```
run_tags.dox_aba_matrix(directory, rr, yy, ry, replicate_number)
```
Additionally, if access to a parallel job creator is available, the main method of run_tags can be run, requiring a job file with each line a specific aba/dox ratio. On the PACE cluster, this is the command pace-gnu-job.

# Questions 
For any questions about the code, contact Andrew Jin at ajin40 [at] gatech.edu
