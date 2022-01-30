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
packaged installation of HOOMD-blue and OpenMM. It is custom and comes with the plugins that are needed for PySAGES for advanced sampling. We copy it from google drive and install pysages for it. We also have a google collab that performs this installation for reference. This is however meant as an explanation of how to install these tools in your environment.

For your work, you can use the same environment if you want to work with google colab or install HOOMD-blue with the installation instructions from the [documentation](https://hoomd-blue.readthedocs.io/en/stable/installation.html) and use local Jupyter notebooks or python scripts.

<!-- #endregion -->

```bash id="nMThqa-DjVcb"

BASE_URL="https://drive.google.com/u/0/uc?id=1hsKkKtdxZTVfHKgqVF6qV2e-4SShmhr7&export=download"
wget -q --load-cookies /tmp/cookies.txt "$BASE_URL&confirm=$(wget -q --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate $BASE_URL -O- | sed -rn 's/.*confirm=(\w+).*/\1\n/p')" -O pysages-env.zip
rm -rf /tmp/cookies.txt
```

```python colab={"base_uri": "https://localhost:8080/"} id="25H3kl03wzJe" outputId="14f0f14d-946d-42d3-ae99-974c3c5e87a9"
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

<!-- #region id="bXVUTurBjFaH" -->
# HOOMD-blue

For this part of the lecture, we introduce common molecular dynamics simulation engines. We start with HOOMD-blue which is an engine centered around the use of GPUs and developed with a python frontend. We are going to use version 2.X for this lecture. There is a new version available 3.X soon, that integrates even tighter with Python.

The homepage of HOOMD-blue can be found [here](https://glotzerlab.engin.umich.edu/hoomd-blue/). It is actively developed on github.com [here](https://github.com/glotzerlab/hoomd-blue) and has an excellent documentation [here](https://hoomd-blue.readthedocs.io/en/stable/package-hoomd.html).
Take some time and start exploring the documentation. It has good examples and explains different parts of the code well. In the context of this course it makes to focus on the core `hoomd` and the `hoomd.md` package.
<!-- #endregion -->

```python id="Qrd0g908yFZt"
import hoomd
import hoomd.md
import numpy as np
```

<!-- #region id="Y0nUxjlkzxFj" -->
HOOMD-blue 2 operates with fixed context for the entire runtime, which also determines the available hardware. So our first step is to initialize this context.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="EVqV_81izwXw" outputId="6d4a162f-8a4b-433a-a897-88c2a973dd6d"
hoomd.context.initialize("")
```

<!-- #region id="b-lE1uxP0DbQ" -->
The output informs us about version information, essential for reproducibility.
And who are the authors of the software and which paper we should cite if we use HOOMD-blue for research.

After that, we are informed which hardware could be detected and HOOMD-blue is running on.
<!-- #endregion -->

<!-- #region id="M4X20ddQxcbZ" -->
## Simulation of soft, coarse-grained polymers

We start with a familiar system. A melt of polymer chains densely fills the simulation box. The first step is to obtain the initial position of all chains.
One way to obtain this for HOOMD-blue is to generate a simulation snapshot and populate the snapshot with the values you are interested in.

Here we set up a system of 25 polymer chains with chain lengths $N=35$ each.

### Generation of initial conditions
<!-- #endregion -->

```python id="6sqwsd6ZzMlI"
n = 25
N = 35
L = 4
b0 = 0.75
snapshot = hoomd.data.make_snapshot(N=n*N,
                                    box=hoomd.data.boxdim(Lx=L, Ly=L, Lz=L),
                                    particle_types=['A'],
                                    bond_types=['backbone'])
```

<!-- #region id="5elC4Pd-0jXx" -->
Now we set up the initial position of each of the polymers. We do that, by randomly selecting the position of the first bead and growing the chain with the ideal gas distribution from there.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="TKWK-eQM0i_A" outputId="12bf38b5-71e0-4474-b39d-8ba35541c281"
np_seed = 42
np_rng = np.random.default_rng(np_seed)

np_pos = np.zeros((n,N, 3))
for poly in range(n):
  np_pos[poly, 0] = np_rng.uniform(-L/2, L/2, (3,))
  for mono in range(1, N):
    # add gaussian distribution to previous bead pos.
    np_pos[poly, mono] = np_pos[poly, mono-1] + np_rng.normal(0, b0, (3,))
# HOOMD-blue does not know about polymers, so we reshape the positions to be a single array
np_pos = np_pos.reshape((n*N, 3))
print(np_pos)
```

<!-- #region id="5G7DW6qx29sj" -->
These positions are a great initial condition for beads, however, we have not respected the periodic boundary conditions with them so far. So we correct the periodic boundary conditions and remember which periodic image each particle belongs to.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="A8UEgKu53U_n" outputId="0ee594fd-d263-4ac4-f535-afef52194474"
box_L = np.asarray([L, L, L])
np_image = np.rint(np_pos / box_L).astype(int)
np_pos -= np_image * box_L

snapshot.particles.position[:] = np_pos
snapshot.particles.image[:] = np_image

print(snapshot.particles.image)
print(snapshot.particles.position)
```

<!-- #region id="ATDHIqRl4jwN" -->
We see that the first polymer was well within the box, but the last polymer was partly outside the box and has been folded back into the box.

Next are the bonds. We need to connect all the beads along with the individual polymers. In total we have $n(N-1)$ bonds.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="5dnvlDRV5alR" outputId="d19b7749-5fa9-48c1-f20d-009e0e405d2b"
snapshot.bonds.resize(n*(N-1))
# All bonds have the same interatin parameter, so they are of type `backbone` = 0
snapshot.bonds.typeid[:] = np.zeros(snapshot.bonds.N, dtype=int)
# Now we add the bonds as necessary
bonds = []
for poly in range(n):
  for mono in range(1, N):
    bond = np.asarray([mono-1, mono], dtype=int)
    # Shift to global particle index
    bond += poly*N
    bonds.append(bond)

snapshot.bonds.group[:] = bonds
print(snapshot.bonds.group)
```

<!-- #region id="ntsNkRRA7YP5" -->
As initial velocity, we choose the Maxwell-Boltzmann distribution.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="IpAICi3H7XqI" outputId="7d297603-a2f4-4238-dde7-0c5f089694ae"
snapshot.particles.velocity[:] = np_rng.normal(0, 1., (n*N, 3))
print(snapshot.particles.velocity)
```

<!-- #region id="hyS8eE2g7vfr" -->
We do not have to assign masses, since we leave all particles with the standard mass of 1.

### Simulation Setup
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="OsEBNyn19ApZ" outputId="69e1e069-191c-4bdf-a048-9e8166dcedae"
system = hoomd.init.read_snapshot(snapshot)
```

<!-- #region id="2kt0pmis-RGQ" -->
For effective calculation of non-bonded forces, we need a Verlet-Neighbor list. We use the standard list here, which uses a grid for fast computation.
<!-- #endregion -->

```python id="NglMGhnm-fRh"
nl = hoomd.md.nlist.cell()
```

<!-- #region id="QQyHVJuR9Ifs" -->
We need to define the interaction between the particles.
As non-bonded forces with choose the soft DPD potential, which has the shape $V(r) = A\cdot (r_\text{cut} - r) - 0.5 A/r_\text{cut} (r_\text{cut}^2 - r^2)$ for distances smaller than $r_\text{cut}$ and zero otherwise.
We use the conservative form here because we don't want the DPD thermostat, instead, we choose Nose-Hover later for the NVT ensemble. 
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="v77UDb8N9HrE" outputId="5b6bfc2b-6c84-4276-9685-ada2dad8a72c"
dpd = hoomd.md.pair.dpd_conservative(r_cut=1, nlist=nl)
dpd.pair_coeff.set("A", "A", A=5.0)
```

<!-- #region id="gt7UGtwr45s_" -->
The next step is to setup a bonded potential for the interaction between beads along the chain. We choose the harmonic potential $V(r) = k/2 (r-r_0)^2$ with a resting length $r_0 = 0$ and k that corresponds to the chosen $k = 3/b_0^2$.
<!-- #endregion -->

```python id="LSgCUDyV45Cs"
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set("backbone", k = 3./b0**2, r0=0)
```

<!-- #region id="NRqeha1P7mSZ" -->
Next, we choose an integrator and the time step. As an integrator, we use NVT integration, which means for HOOMD-blue two-step Velocity-Verlet integration with a Nose-Hover thermostat.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="gEfD8s0H88b7" outputId="07b2ccea-0d6c-4479-9c5c-f317f91395fe"
dt=1e-3
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.nvt(hoomd.group.all(), kT=1.0, tau=1.0)
```

<!-- #region id="-135fr8q9gNT" -->
We do not only want to run the simulation, but we also would like to measure some properties during the simulation run. Hoomd has a logging system implemented that allows you to write observables to disk during the simulation run. Hoomd-blue has common observables implemented but lets user also specify their own observables. 
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="-QWo5uDh-k3O" outputId="9321ea5a-e68e-43c4-eac6-360c034654ca"
quantities_request = ['time', 'temperature',"pressure",'potential_energy','kinetic_energy', "bond_harmonic_energy"]
hoomd.analyze.log("obs.txt", quantities_request, period=100, overwrite=True)
```

<!-- #region id="JASGrFHD9gKJ" -->
We also want to log a trajectory, that contains all particle positions, bond information, etc. This is primarily useful for visualization and post-processing for more observables. We write the data in Hoomd's data format [gsd](https://gsd.readthedocs.io/en/stable/).
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="S6NRmrS0_mWO" outputId="738b2032-3ae9-41ce-d58f-4abb97b4f9a0"
hoomd.dump.gsd("traj.gsd", 1e3, hoomd.group.all(), overwrite=True)
```

<!-- #region id="PuZhE13gALGj" -->
Now we can simulate a given number of time steps.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="jdFVKJhAAKnx" outputId="83c0ef09-096f-4e75-ec97-4e564c187ba1"
timesteps = 5e5
hoomd.run(timesteps)
```

<!-- #region id="vItDwIP8DPdw" -->
At the end of the run, it is best practice to store the last configuration. This allows you to continue the simulation from the last step if it turns out, that you require more data for your analysis. Or if your wall time limitations on supercomputer resources.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="3Y82L0E0SJcb" outputId="75e5641c-ac92-4e8d-b543-72dc15ba3c5b"
hoomd.dump.gsd("final.gsd", None, hoomd.group.all(), overwrite=True, truncate=True)
```

<!-- #region id="jaiNo0WJToYH" -->
# Analysis

## Visualization

We can use tools like Ovito[link text](https://www.ovito.org/) to visualize and inspect the full trajectory. Ovito has a free, open-source version that allows basic functions. A pro version is available but covers only functions that are beyond the scope of this course.

Ovito does not integrate seamlessly into this notebook, so we doing this offline.

## Plotting of data

However, during the simulation run, we stored some quantities in plain text format and we can use numpy and matplotlib to visualize this data.

<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="BwnmYH-8Uyt8" outputId="5408a17e-a3b0-4d21-e8ed-01869cdf4686"
!head obs.txt
import matplotlib.pyplot as plt

timestep, time, temp, pressure, Epot, Ekin, Ebond = np.loadtxt("obs.txt", skiprows=1, unpack=True)
```

<!-- #region id="c9cI7HqIXJBJ" -->
### Temperature

The Nose-Hover thermostat is supposed to keep the temperature at $k_BT=1$, so we can see if that is the case.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 301} id="pbLOtRfwW0Bt" outputId="d93e7879-1e78-4834-8b93-f663d098e2c4"
fig, ax = plt.subplots()
ax.set_xlabel(r"time [$\tau$]")
ax.set_ylabel(r"$T$ [$\epsilon$]")
ax.plot(time, temp, label="instantanenous temeprature")
```

<!-- #region id="hjxu9uRHX2lF" -->
The observed oscillations are very common for the Nose-Hover algorithm. The time scale of these oscillations should be on the time scale of our friction time of the thermostat $\tau_\text{NH} =1 \tau$. And it is very crucial to equilibrate the system and its temperature before any measurements. Even as the oscillations do not decay completely in this example, we see, that the average temperature we set for the thermostat is achieved.

### Pressure

Since we perform NVT simulations the thermodynamic variable to volume $V$ pressure $p$ is a result of the simulations.
We can visualize this as a function of time.

<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 301} id="2IA0gb-CXy00" outputId="114a5fda-6f0a-4430-8493-8914cd02fa7e"
fig, ax = plt.subplots()
ax.set_xlabel(r"time [$\tau$]")
ax.set_ylabel(r"$p$ [$\epsilon/\sigma^3$]")
ax.plot(time, pressure, label="pressure")
print(pressure.mean(), pressure.std())
```

<!-- #region id="rnCOga5maegP" -->
We also observe the unphysical effects of equilibration here at first. After that, we see how the pressure converges to a mean value. We expect this to be a thermal average of the pressure $\langle p \rangle$.

I print here simply the mean and std for reference for a real analysis we need to remove the equilibration period first and for the error estimation we need to perform a block analysis.

### Energy analysis

The last analysis we are interested in here is the evolution of energy.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 301} id="S2TtsDP0aa-x" outputId="7df548db-2849-47e2-ab63-3d63ac6e9ae7"
fig, ax = plt.subplots()
ax.set_xlabel(r"time [$\tau$]")
ax.set_ylabel(r"$E$ [$\epsilon$]")
ax.plot(time, Epot, label=r"$E_\mathrm{pot.}$")
ax.plot(time, Ekin, label=r"$E_\mathrm{kin.}$")
ax.plot(time, Ebond, label=r"$E_\mathrm{bond}$")
ax.legend(loc="best")
```

<!-- #region id="Vo-o5ns-cAqC" -->
We see that we start with a very high potential energy, also driven by overstretched bonds. This gets naturally traded for kinetic energy as the forces accelerate the particles. Here kicks the thermostat in and sheds the kinetic energy. This happens in oscillation until the equilibrium average potential and kinetic energies are reached.
<!-- #endregion -->
