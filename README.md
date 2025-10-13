# BTBD (Block-Tether Brownian Dynamics)

Inputs: design file, which specifies the geometry of the rigid bodies, the construction of origamis from these rigid bodies, and the linkers that connect origamis.

Output: geometry and input files for running the simulation in LAMMPS.

See the main script (BTBD.m) for a detailed description of the code, including how to construct design files and adjust simulation parameters.

## Input File

Blocks: rigid bodies
- name - how to identify the block for patches and origamis
- pattern - arangement of beads in the xy-plane (see block object)
- height - number of beads along the z-axis

Patches: locations on blocks used for bonded interactions
- name - how to identify the patch for connections and linkers
- block - name of block on which to place the patch
- x,y,z - patch location in block coordinate system

Origami: collection of connected blocks
- name - how to identify the origami for defining parameters and adding linkers
- count - number of origamis to create
- blocks - names of block types (blocks are indexed in the given order)
- conn - permanent bond between between indexed blocks at named patches

Linker: breakable bond
- origami - name of origamis on which to create the linker
- block indices - index, or "A" for all blocks, or "B" for all but the last block
- location - where to place the ends of the linker