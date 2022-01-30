---
jupyter:
  jupytext:
    formats: ipynb,md
    main_language: python
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.6
  kernelspec:
    display_name: Python 3
    name: python3
---

<!-- #region id="p49wJ0IjLAVD" -->

# Setup of the environment
<!-- #endregion -->

<!-- #region id="WM_9PpDwKuoA" -->
First, we are setting up our environment. We use an already compiled and 
packaged installation of HOOMD-blue and OpenMM. It is custom and comes with the plugins that are needed for PySAGES for advanced sampling. We copy it from google drive and install pysages for it. We also have a google collab that performs this installation for reference. This is however meant as an explanation of how to install these tools in your own environment.

For your own work, you can use the same environment if you want to work with google colab or install HOOMD-blue with the installation instructions from the [documentation](https://hoomd-blue.readthedocs.io/en/stable/installation.html) and use local Jupyter notebooks or python scripts.

<!-- #endregion -->

```bash id="nMThqa-DjVcb"

BASE_URL="https://drive.google.com/u/0/uc?id=1hsKkKtdxZTVfHKgqVF6qV2e-4SShmhr7&export=download"
wget -q --load-cookies /tmp/cookies.txt "$BASE_URL&confirm=$(wget -q --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate $BASE_URL -O- | sed -rn 's/.*confirm=(\w+).*/\1\n/p')" -O pysages-env.zip
rm -rf /tmp/cookies.txt
```

```python colab={"base_uri": "https://localhost:8080/"} id="25H3kl03wzJe" outputId="e0f68c85-9253-43b5-f68a-392393e86e13"
%env PYSAGES_ENV=/env/pysages
```

```bash id="CPkgxfj6w4te"

mkdir -p $PYSAGES_ENV
unzip -qquo pysages-env.zip -d $PYSAGES_ENV
```

```python id="JMO5fiRTxAWB"
import os
import sys

ver = sys.version_info

sys.path.append(os.environ["PYSAGES_ENV"] + "/lib/python" + str(ver.major) + "." + str(ver.minor) + "/site-packages/")
```

```python id="Qrd0g908yFZt" colab={"base_uri": "https://localhost:8080/"} outputId="de404695-1c7f-4327-c301-03204a43653f"
import hoomd
import hoomd.md
import numpy as np
from numpy import sqrt, pi
hoomd.context.initialize("")
```

<!-- #region id="Y0nUxjlkzxFj" -->
HOOMD-blue 2 operates with fixed context for the entire runtime, which also determines the available hardware. So our first step is to initialize this context.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="EVqV_81izwXw" outputId="11542cd7-b51f-4931-c9b4-9be96dc04378"
hoomd.context.initialize("")
```

<!-- #region id="b-lE1uxP0DbQ" -->
The output informs us about version information, essential for reproducibility.
And who are the authors of the software and which paper we should cite if we useHOOMD-blue for research.

After that, we are informed which hardware could be detected and HOOMD-blue is running on.
<!-- #endregion -->

<!-- #region id="v97i6xNLno0D" -->
# Butane simulation

After we explored a system of abstract/coarse-grained nature we are investigating now an all-atom system of [Butane](https://en.wikipedia.org/wiki/Butane).
This tutorial is inspired by the original [SSAGES tutorial](https://github.com/SSAGESproject/SSAGES/tree/release-0.9/Examples/User/Umbrella/HOOMD) and adjust the SSAGES tutorial for PySAGES.

The biggest difference between such a detailed system is that the force-field of all-atom simulations is more complex. Each bead needs to be identified with the appropriate force-field parameter and the interactions (especially non-bonded interactions) are stiffer than before. This requires more care when setting up the initial conditions to avoid overlaps of the individual atoms.

This tutorial is meant as a demonstration. Do not use this forcefield for research. We will later see how accurate forcefields can be obtained.

Here we set up the initial conditions from a HOOMD-blue snapshot.
<!-- #endregion -->

```python id="pQA3DdRepwGJ"
snapshot = hoomd.data.make_snapshot(N=14,
                                    box=hoomd.data.boxdim(Lx=41, Ly=41, Lz=41),
                                    particle_types=['C', 'H'],
                                    bond_types=['CC', 'CH'],
                                    angle_types=['CCC', 'CCH', 'HCH'],
                                    dihedral_types=['CCCC', 'HCCC', 'HCCH'],
                                    pair_types=['CCCC', 'HCCC', 'HCCH'])
```

<!-- #region id="BccxD_0w_3EP" -->
This snapshot contains 14 particles and we have to assign which particle represents a C-atom and which particle is an H-atom.
<!-- #endregion -->

```python id="EVU5jqoG_4pR"
snapshot.particles.typeid[0] = 0
snapshot.particles.typeid[1:4] = 1
snapshot.particles.typeid[4] = 0
snapshot.particles.typeid[5:7] = 1
snapshot.particles.typeid[7] = 0
snapshot.particles.typeid[8:10] = 1
snapshot.particles.typeid[10] = 0
snapshot.particles.typeid[11:14] = 1


```

<!-- #region id="HeJNkA4GAcgy" -->
We also assign each of the atoms an initial condition. Information about the good initial condition for molecules can be generated for example via [RDkit](https://rdkit.org/) or for a biological system from the [protein databank](https://www.rcsb.org/), or via a small custom MC simulation. Here we are fortunate to have some valid initial positions.
<!-- #endregion -->

```python id="rAdVHbdtAXeg"
  positions = np.array([
      [-2.990196,  0.097881,  0.000091],
      [-2.634894, -0.911406,  0.001002],
      [-2.632173,  0.601251, -0.873601],
      [-4.060195,  0.099327, -0.000736],
      [-2.476854,  0.823942,  1.257436],
      [-2.832157,  1.833228,  1.256526],
      [-2.834877,  0.320572,  2.131128],
      [-0.936856,  0.821861,  1.258628],
      [-0.578833,  1.325231,  0.384935],
      [-0.581553, -0.187426,  1.259538],
      [-0.423514,  1.547922,  2.515972],
      [-0.781537,  1.044552,  3.389664],
      [ 0.646485,  1.546476,  2.516800],
      [-0.778816,  2.557208,  2.515062]])

  reference_box_low_coords = np.array([-22.206855, -19.677099, -19.241968])
  box_low_coords = np.array([-snapshot.box.Lx/2,
                             -snapshot.box.Ly/2,
                             -snapshot.box.Lz/2])
  positions += (box_low_coords - reference_box_low_coords)

  snapshot.particles.position[:] = positions[:]
  
```

<!-- #region id="i4m0sVz8AXUk" -->
Next are masses and charges. Note how this automatically chooses the unit system for the simulations.
<!-- #endregion -->

```python id="kmhU4IkZNBl4"
mC = 12.00
mH = 1.008
snapshot.particles.mass[:] = [mC, mH, mH, mH,
                              mC, mH, mH,
                              mC, mH, mH,
                              mC, mH, mH, mH]
reference_charges = np.array([-0.180000, 0.060000, 0.060000, 0.060000,
                              -0.120000, 0.060000, 0.060000,
                              -0.120000, 0.060000, 0.060000,
                              -0.180000, 0.060000, 0.060000, 0.060000])
charge_conversion = 18.22262
snapshot.particles.charge[:] = charge_conversion * reference_charges[:]
```

<!-- #region id="IvRvQXLqq6jt" -->
The next step is to set up the bond potentials and assign the different bond types to the individual bonds. We are also using angle potentials, dihedrals, and special pairs for the interaction of this molecule. See the hoomd documentation for details.
<!-- #endregion -->

```python id="fgqge-dKM9bB"
snapshot.bonds.resize(13)
snapshot.bonds.typeid[0:3] = 1
snapshot.bonds.typeid[3] = 0
snapshot.bonds.typeid[4:6] = 1
snapshot.bonds.typeid[6] = 0
snapshot.bonds.typeid[7:9] = 1
snapshot.bonds.typeid[9] = 0
snapshot.bonds.typeid[10:13] = 1

snapshot.bonds.group[:] = [[0, 2], [0, 1], [0, 3], [0, 4],
                           [4, 5], [4, 6], [4, 7],
                           [7, 8], [7, 9], [7, 10],
                           [10, 11], [10, 12], [10, 13]]

snapshot.angles.resize(24)
snapshot.angles.typeid[0:2] = 2
snapshot.angles.typeid[2] = 1
snapshot.angles.typeid[3] = 2
snapshot.angles.typeid[4:8] = 1
snapshot.angles.typeid[8] = 0
snapshot.angles.typeid[9] = 2
snapshot.angles.typeid[10:14] = 1
snapshot.angles.typeid[14] = 0
snapshot.angles.typeid[15] = 2
snapshot.angles.typeid[16:21] = 1
snapshot.angles.typeid[21:24] = 2

snapshot.angles.group[:] = [[1, 0, 2], [2, 0, 3], [2, 0, 4],
                            [1, 0, 3], [1, 0, 4], [3, 0, 4],
                            [0, 4, 5], [0, 4, 6], [0, 4, 7],
                            [5, 4, 6], [5, 4, 7], [6, 4, 7],
                            [4, 7, 8], [4, 7, 9], [4, 7, 10],
                            [8, 7, 9], [8, 7, 10], [9, 7, 10],
                            [7, 10, 11], [7, 10, 12], [7, 10, 13],
                            [11, 10, 12], [11, 10, 13], [12, 10, 13]]

snapshot.dihedrals.resize(27)
snapshot.dihedrals.typeid[0:2] = 2
snapshot.dihedrals.typeid[2] = 1
snapshot.dihedrals.typeid[3:5] = 2
snapshot.dihedrals.typeid[5] = 1
snapshot.dihedrals.typeid[6:8] = 2
snapshot.dihedrals.typeid[8:11] = 1
snapshot.dihedrals.typeid[11] = 0
snapshot.dihedrals.typeid[12:14] = 2
snapshot.dihedrals.typeid[14] = 1
snapshot.dihedrals.typeid[15:17] = 2
snapshot.dihedrals.typeid[17:21] = 1
snapshot.dihedrals.typeid[21:27] = 2

snapshot.dihedrals.group[:] = [[2, 0, 4, 5], [2, 0, 4, 6], [2, 0, 4, 7],
                               [1, 0, 4, 5], [1, 0, 4, 6], [1, 0, 4, 7],
                               [3, 0, 4, 5], [3, 0, 4, 6], [3, 0, 4, 7],
                               [0, 4, 7, 8], [0, 4, 7, 9], [0, 4, 7, 10],
                               [5, 4, 7, 8], [5, 4, 7, 9], [5, 4, 7, 10],
                               [6, 4, 7, 8], [6, 4, 7, 9], [6, 4, 7, 10],
                               [4, 7, 10, 11], [4, 7, 10, 12], [4, 7, 10, 13],
                               [8, 7, 10, 11], [8, 7, 10, 12], [8, 7, 10, 13],
                               [9, 7, 10, 11], [9, 7, 10, 12], [9, 7, 10, 13]]

snapshot.pairs.resize(27)
snapshot.pairs.typeid[0:1] = 0
snapshot.pairs.typeid[1:11] = 1
snapshot.pairs.typeid[11:27] = 2
snapshot.pairs.group[:] = [
    # CCCC
    [0, 10],
    # HCCC
    [0, 8], [0, 9], [5, 10], [6, 10],
    [1, 7], [2, 7], [3, 7],
    [11, 4], [12, 4], [13, 4],
    # HCCH
    [1, 5], [1, 6], [2, 5], [2, 6], [3, 5], [3, 6],
    [5, 8], [6, 8], [5, 9], [6, 9],
    [8, 11], [8, 12], [8, 13], [9, 11], [9, 12], [9, 13]
]
```

<!-- #region id="W9nDKFWrsPJL" -->
## Simulation setup

As before after we obtained the snapshot we are initializing the actual simulation setup.
Because of the more complicated nature of all-atom butane, we need a more complex forcefield. The potentials and their parameter are part of the forcefield. For now, we just accept them, we will later see how to obtain actual forcefields for molecules.
<!-- #endregion -->

```python id="onu9MNuXzQNT" colab={"base_uri": "https://localhost:8080/"} outputId="e9aa966e-dbac-40d7-9c02-5b79d2e336dd"
system = hoomd.init.read_snapshot(snapshot)

#Neighbor list, reset some of HOOMD's automatic exclusion rules.
nl_ex = hoomd.md.nlist.cell()
nl_ex.reset_exclusions(exclusions=['1-2', '1-3', '1-4'])

# Base non-bonded interactions: Lennard-Jones, that parameters are in simulation units and part of the forcefield.
lj = hoomd.md.pair.lj(r_cut=12.0, nlist=nl_ex)
lj.pair_coeff.set('C', 'C', epsilon=0.07, sigma=3.55)
lj.pair_coeff.set('H', 'H', epsilon=0.03, sigma=2.42)
lj.pair_coeff.set('C', 'H', epsilon=sqrt(0.07*0.03), sigma=sqrt(3.55*2.42))

# For all-atom system the partial charge interactions are important.
coulomb = hoomd.md.charge.pppm(hoomd.group.charged(), nlist=nl_ex)
# Parameters for the Ewald summation
coulomb.set_params(Nx=64, Ny=64, Nz=64, order=6, rcut=12.0)

# Bond potentials are harmonic
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('CC', k=2*268.0, r0=1.529)
harmonic.bond_coeff.set('CH', k=2*340.0, r0=1.09)

# Angle potentials
angle = hoomd.md.angle.harmonic()
angle.angle_coeff.set('CCC', k=2*58.35, t0=112.7 * pi / 180)
angle.angle_coeff.set('CCH', k=2*37.5, t0=110.7 * pi / 180)
angle.angle_coeff.set('HCH', k=2*33.0, t0=107.8 * pi / 180)

# Dihedrals 
dihedral = hoomd.md.dihedral.opls()
dihedral.dihedral_coeff.set('CCCC', k1=1.3, k2=-0.05, k3=0.2, k4=0.0)
dihedral.dihedral_coeff.set('HCCC', k1=0.0, k2=0.0, k3=0.3, k4=0.0)
dihedral.dihedral_coeff.set('HCCH', k1=0.0, k2=0.0, k3=0.3, k4=0.0)

# These are special LJ particle pairs, that are treated like bonded. LJ just better performance
lj_special_pairs = hoomd.md.special_pair.lj()
lj_special_pairs.pair_coeff.set('CCCC', epsilon=0.07, sigma=3.55, r_cut=12.0)
lj_special_pairs.pair_coeff.set('HCCH', epsilon=0.03, sigma=2.42, r_cut=12.0)
lj_special_pairs.pair_coeff.set('HCCC', epsilon=sqrt(0.07 * 0.03),
                                sigma=sqrt(3.55*2.42), r_cut=12.0)

coulomb_special_pairs = hoomd.md.special_pair.coulomb()
coulomb_special_pairs.pair_coeff.set('CCCC', alpha=0.5, r_cut=12.0)
coulomb_special_pairs.pair_coeff.set('HCCC', alpha=0.5, r_cut=12.0)
coulomb_special_pairs.pair_coeff.set('HCCH', alpha=0.5, r_cut=12.0)
```

<!-- #region id="E04RypAKIa_X" -->
Next is the integrator. We choose Nose-Hover again.
<!-- #endregion -->

```python id="AakEHD9lJ2R_"
kT = 0.596161
dt = 0.02045
hoomd.md.integrate.mode_standard(dt=dt)
integrator = hoomd.md.integrate.nvt(group=hoomd.group.all(), kT=kT, tau=100*dt)
```

<!-- #region id="QNsLGYvaJ8jA" -->
We also would like to log some quantities calculated on the fly and the trajectory.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="s0DVXEX_J43N" outputId="0e0136b2-6da0-4eb8-fc04-59630e247188"
quantities_request = ['time', 'temperature',"pressure",'potential_energy','kinetic_energy']
hoomd.analyze.log("obs.txt", quantities_request, period=100, overwrite=True)
hoomd.dump.gsd("traj.gsd", 1e3, hoomd.group.all(), overwrite=True)
```

<!-- #region id="m8cRkXV_K5F-" -->
And now, we run the simultions.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="EXlGGmk5K11k" outputId="7464844e-81fa-4ec6-a4b1-d76bce34924a"
timesteps = 1e5
hoomd.run(timesteps)
```

<!-- #region id="y_FMMdaqODOS" -->
And from here we can visualize and analyze the trajectory.

However, setting the forcefields up manually is a challenging task that we address in more detail next.
<!-- #endregion -->
