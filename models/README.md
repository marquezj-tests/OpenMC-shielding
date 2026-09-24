# Models

## `shultis_deep_shielding.py`

OpenMC reproduction of the deep-shielding example in Shultis & Faw,
*An MCNP Primer* (2011), Section 4.4.3 "An Example" / Figure 5 (p. 31) —
the book's illustration of why deep-penetration shielding problems need
variance reduction (geometry splitting or weight windows).

**Problem:** a point isotropic 7-MeV photon source at the center of a
30-cm-thick natural-iron spherical shell (r = 30-60 cm, 7.86 g/cm³). The
ambient dose equivalent H*(10) is tallied on a sphere at r = 160 cm, using
the ICRP-1987 fluence-to-dose conversion factors given on the book's
DE2/DF2 cards. The book states the unshielded (no shield) dose is
6.013e-17 Sv/photon.

**Run it:**

```bash
export OPENMC_CROSS_SECTIONS=~/nuclear_data/cross_sections.xml

python models/shultis_deep_shielding.py              # shielded case
python models/shultis_deep_shielding.py --no-shield   # bare-source validation
```

**Validation:** the `--no-shield` case reproduces the book's own quoted
reference dose (6.013e-17 Sv/photon) to within 0.2% (5.999e-17 Sv/photon
here), confirming the geometry, source, and dose-response tally are
correctly implemented independent of the iron cross-section data.

**Shielded result:** with 1,000,000 analog histories (20 batches x 50,000
particles, no variance reduction), the model gives an ambient dose of
~1.78e-19 +/- 0.03e-19 Sv/photon at r = 160 cm — about a factor of 337
attenuation from the bare-source case, consistent with roughly 6 mean
free paths of Compton-dominated attenuation in iron at 7 MeV plus
buildup. This converges comfortably on modern hardware; the book's own
literal `nps 10000` (10,000 histories) would give much poorer statistics,
which is exactly why it uses this problem to motivate geometry splitting
and weight windows as the next topic.

**Tally technique:** MCNP's F2 (surface flux) tally has no direct OpenMC
equivalent, so the model uses the standard substitute: a thin (1-cm)
spherical shell at r = 160 cm with a track-length `flux` score divided by
the shell volume, combined with an `EnergyFunctionFilter` built from the
DE2/DF2 table (linear-linear interpolation, matching MCNP's default) to
reproduce the dose-weighted tally.

**Known deviations from the MCNP physics options** (`phys:p 100 1 1`):
- *No bremsstrahlung*: reproduced via `settings.electron_treatment = 'led'`.
- *No coherent (Rayleigh) scattering*: OpenMC has no per-run toggle for
  this; Rayleigh scattering is included whenever the photon data provides
  it. The effect on this problem is small (forward-peaked, doesn't
  materially change deep penetration).

The shield is pre-divided into the ten 3-cm sub-shells the book uses for
its geometry-splitting scheme (cell importances 1, 2, 4, ..., 512), so
weight windows or another variance-reduction layer can be added later
without changing the geometry. This baseline model itself runs analog
(no biasing).
