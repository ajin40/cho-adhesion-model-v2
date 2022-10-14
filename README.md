# Vicsek-like ABM of cell self-organization in cadherin-linked synthetic switch

Link to manuscript: To be added

## Title
Synthetic symmetry breaking and programmable multicellular structure formation
## Authors
Noreen Wauford, Akshay Patel, Jesse Tordoff, Casper Enghuus, Andrew Jin<sup>+</sup>, Jack Toppen, Melissa L. Kemp, Ron Weiss

<sup>+</sup>: Code author

# Model Description
Cells are modeled as rigid spheres of three cell types: blue, red, and yellow cells all with equal cell diameter σ. The three cell types vary in adhesion strength, and are in ascending order: blue, red, and yellow, with yellow-yellow adhesion the highest. All cells are initialized as blue cells with a fixed probability of transitioning to a cadherin expressing (red, yellow) state determined by the final fluorescence intensity expression of the corresponding Dox/ABA dosage experiment, p_total cells committed. The fixed transition probability at each time step (x) is back calculated from the final fluorescence intensity expression:

$$(1-x)^72h=(1-p_total cells committed)$$

When a cell is indicated to transition to a cadherin expressing state, the probability of choosing the red or yellow states is also determined by the fluorescence intensity expression. The probability of choosing the one of the states is equal to the experimental intensity expression ratio of the respective cells over the total non-blue intensity expression. 

Cell locations are stored in a NumPy array and used for local cell interaction calculations. Cell locations are updated each time step (10 simulated seconds), at a fixed velocity of 0.3 * σ (cell radius). The direction of each cell’s velocity vector is determined by the normalized sum of all forces acting upon each cell. These forces are Brownian motion, a gravity field, and the local cell-cell interactions of homotypic and heterotypic adhesion and steric repulsion. 

All cells experience Brownian motion and a gravity field designed to mimic cell behavior in an ultra-low binding U-shaped bottom well. To represent the curvature of the well, the magnitude of the gravity field is dependent on the xy planar distance $$(|r_(i,xy)|)$$ from the center of the simulation space. U_g represents the strength of the force and R as the radius of the well.

$$f_i=-U_g r_i  (|r_(i,xy)|)/R$$

The presence of the pairwise forces is dependent on the distance $$r_ij$$ between two cells. At $$|r_ij<1|$$, two cells are overlapping, and each cell experiences a strong $$(U_hc=10^4)$$ steric repulsion force to remove overlap of the rigid bodies.

$$f_ij=-U_hc  r_ij/|r_ij|$$ 

At $$|1< r_ij<1.6*σ|$$, cells experience attraction due to cell-cell adhesion. This is described by

$$f_ij=U_ij (r_ij-r_e )  r_ij/|r_ij|$$

Where $$U_ij$$ is the pairwise cell type specific adhesion force. There are 6 cell-type specific adhesion parameters, of which 3 are for the homotypic ($$U_bb,U_rr,U_yy$$) adhesion forces, and 3 are for the heterotypic ($$U_br,U_by,U_ry$$) adhesion forces. The equilibrium distance of intercellular interaction is denoted as $$r_e$$. The isotropic noise of the attraction force is represented with a Gaussian noise vector $$ξ_ij$$ with magnitude α. At  $$|r_ij>1.6*σ|$$, the cells are sufficiently far and do not interact. As cells sufficiently distant from each other do not influence each other, local interactions are identified through a fixed radius search of 1.6 * σ for nearby cell bodies and then subsequently divided into adhesion and repulsion forces. 


