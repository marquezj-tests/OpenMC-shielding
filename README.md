# OpenMC Shielding

This repository is set up to run [OpenMC](https://openmc.org) photon transport
simulations for shielding analysis.

## Installation

OpenMC does not ship on PyPI; it is distributed via conda-forge. Install
[Miniforge](https://github.com/conda-forge/miniforge) (or any conda/mamba
distribution) and create the environment:

```bash
conda env create -f environment.yml
conda activate openmc-env
```

This installs OpenMC 0.16.0 (CPU, no MPI, no DAGMC) along with NJOY2016,
which OpenMC uses for nuclear data processing.

## Nuclear data

Photon transport in OpenMC needs two kinds of ENDF-derived HDF5 data:

- **Incident photon (photoatomic) data** — one file per element, used for
  the actual photon physics (photoelectric, Compton, pair production,
  coherent scattering).
- **Incident neutron data** — one file per isotope. OpenMC requires this
  for *every* material even in photon-only runs, since it resolves each
  nuclide's atomic weight ratio from the neutron sublibrary.

Fetch both and build the `cross_sections.xml` index with:

```bash
./scripts/setup_nuclear_data.sh          # installs to ~/nuclear_data
export OPENMC_CROSS_SECTIONS=~/nuclear_data/cross_sections.xml
```

The script pulls a pre-converted HDF5 mirror of ENDF/B-VIII.0 (NNDC) from
[openmc-data-storage/ENDF-B-VIII.0-NNDC](https://github.com/openmc-data-storage/ENDF-B-VIII.0-NNDC)
on GitHub. This mirror is used instead of the official `nndc.bnl.gov` /
`openmc.org` hosts because those are not reachable from every network
environment. Set `OPENMC_CROSS_SECTIONS` in your shell profile so it
persists across sessions.

Two isotopes (U235, U238) are stored via Git LFS in that mirror and are not
fetched by an anonymous clone; the setup script skips them and prints a
warning. Everything else — 554 isotopes and all 100 elements' photon data —
is downloaded as regular files.

## Verifying the setup

```bash
python scripts/verify_photon_transport.py
```

This runs a small fixed-source photon transport calculation (1 MeV photons
into a 10 cm lead sphere) and prints the resulting flux and heating tallies.

## Models

- [`models/shultis_deep_shielding.py`](models/shultis_deep_shielding.py) —
  the deep-shielding (point 7-MeV photon source in a 30-cm iron shell)
  example from Shultis & Faw, *An MCNP Primer*. See
  [`models/README.md`](models/README.md).
