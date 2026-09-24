"""OpenMC model of the deep-shielding example from Shultis & Faw,
"An MCNP Primer" (2011), Section 4.4.3 "An Example" / Figure 5, p. 31.

Problem statement (as given in the book): a point isotropic source emitting
7-MeV photons sits at the center of a 30-cm-thick natural-iron spherical
shell (inner radius 30 cm, outer radius 60 cm, density 7.86 g/cm^3). The
ambient dose equivalent H*(10) is sought on a sphere 160 cm from the source
(i.e. 100 cm beyond the outer shield surface). Without the shield, the book
states the dose is 6.013e-17 Sv per source photon.

This module reproduces the MCNP "analog base case" input deck one-for-one:

    10  0        -10           imp:p=1   $ inside of shield
    20  1 -7.86  10 -20        imp:p=1   $ iron shell
    30  0        20 -50        imp:p=1   $ void outside shld, inside detect
    40  0        50 -100       imp:p=1   $ void past detector
    50  0        100           imp:p=0   $ vacuum outside problem boundary

    10 so 30.0        $ inner shield surface
    20 so 60.0         $ outer shield surface
    50 so 160.0        $ detector surface
    100 so 10.E+02      $ spherical problem boundary (10 m)

    SDEF erg=7.00 par=2      $ 7-MeV point isotropic photon source
    phys:p 100 1 1           $ no bremsstrahlung; no coherent scattering
    f2:p 50                  $ ambient dose on surface 50
    de2/df2                  $ ICRP-1987 H*(10) fluence-to-dose factors
    m1  26000 -1.00000       $ natural iron, density 7.86 g/cm^3

The 30-cm shield is additionally pre-divided into ten 3-cm sub-shells
(cells 20-29 in the book), which is what the book's geometry-splitting
variance-reduction scheme (cell importances 1,2,4,...,512) is later built
on. The physical model is identical whether or not the shield is
subdivided; the subdivision is kept here so weight windows or cell-based
variance reduction can be layered on later without changing the geometry.

Two deviations from the exact MCNP physics options, both noted because
OpenMC does not expose the same switches:
  - "no coherent scattering" (phys:p ...1): OpenMC includes coherent
    (Rayleigh) scattering whenever the photon data provides it, and there
    is no per-run toggle to disable it. The effect on this problem is
    small (Rayleigh scattering is forward-peaked and does not appreciably
    change energy deposition or deep penetration).
  - "no bremsstrahlung" (phys:p ...1 for electrons): reproduced by using
    settings.electron_treatment = 'led' (local energy deposition), so
    secondary electrons deposit their energy locally instead of producing
    bremsstrahlung photons, matching the book's choice.

The F2 surface-flux tally is reproduced with a thin (1-cm) spherical shell
centered on the detector surface (r=160 cm), a track-length 'flux' score
divided by the shell volume, and an EnergyFunctionFilter built from the
book's DE2/DF2 table (ICRP-1987 ambient dose conversion factors, linear-
linear interpolation, exactly as MCNP's default DE/DF interpolation) to
reproduce the dose-weighted tally.
"""
import argparse

import numpy as np
import openmc

# ICRP-1987 ambient dose equivalent H*(10) conversion factors, exactly as
# given on the DE2/DF2 cards in Figure 5 of the primer.
DE2_MEV = np.array([
    0.100E-01, 0.150E-01, 0.200E-01, 0.300E-01, 0.400E-01, 0.500E-01,
    0.600E-01, 0.800E-01, 0.100E+00, 0.150E+00, 0.200E+00, 0.300E+00,
    0.400E+00, 0.500E+00, 0.600E+00, 0.800E+00, 0.100E+01, 0.150E+01,
    0.200E+01, 0.300E+01, 0.400E+01, 0.500E+01, 0.600E+01, 0.800E+01,
    0.100E+02,
])
DF2_SV_CM2 = np.array([
    0.769E-13, 0.846E-12, 0.101E-11, 0.785E-12, 0.614E-12, 0.526E-12,
    0.504E-12, 0.532E-12, 0.611E-12, 0.890E-12, 0.118E-11, 0.181E-11,
    0.238E-11, 0.289E-11, 0.338E-11, 0.429E-11, 0.511E-11, 0.692E-11,
    0.848E-11, 0.111E-10, 0.133E-10, 0.154E-10, 0.174E-10, 0.212E-10,
    0.252E-10,
])

SOURCE_ENERGY_EV = 7.0e6
SHIELD_INNER_R = 30.0     # cm, surface 10
SHIELD_OUTER_R = 60.0     # cm, surface 20
N_SHIELD_SHELLS = 10      # cells 20-29 in the book
DETECTOR_R = 160.0        # cm, surface 50
DETECTOR_SHELL_HALF_THICKNESS = 0.5   # cm, thin shell around surface 50
OUTER_BOUNDARY_R = 1000.0  # cm (10 m), surface 100
IRON_DENSITY = 7.86        # g/cm^3


def build_model(include_shield: bool = True) -> openmc.Model:
    """Build the OpenMC equivalent of Figure 5.

    include_shield=False reproduces the book's "without shield" reference
    case (dose = 6.013e-17 Sv/gamma) by making the shield cells void, used
    to validate the source, geometry, and dose-tally setup independently
    of the iron cross-section data.
    """
    model = openmc.Model()

    # --- Materials --------------------------------------------------
    iron = openmc.Material(name="natural iron shield")
    iron.add_element("Fe", 1.0)
    iron.set_density("g/cm3", IRON_DENSITY)
    model.materials = openmc.Materials([iron])

    # --- Geometry -----------------------------------------------------
    # Surface 10 / 20 equivalents: ten 3-cm sub-shells spanning 30-60 cm,
    # matching cells 20-29 of the book's geometry-splitting scheme.
    shell_radii = np.linspace(SHIELD_INNER_R, SHIELD_OUTER_R, N_SHIELD_SHELLS + 1)
    shell_surfaces = [openmc.Sphere(r=r) for r in shell_radii]

    cavity = openmc.Cell(
        name="10: inside of shield",
        region=-shell_surfaces[0],
    )

    shield_cells = []
    for i in range(N_SHIELD_SHELLS):
        cell = openmc.Cell(
            name=f"{20 + i}: iron shell {i}",
            region=+shell_surfaces[i] & -shell_surfaces[i + 1],
            fill=iron if include_shield else None,
        )
        shield_cells.append(cell)

    detector_inner = openmc.Sphere(r=DETECTOR_R - DETECTOR_SHELL_HALF_THICKNESS)
    detector_outer = openmc.Sphere(r=DETECTOR_R + DETECTOR_SHELL_HALF_THICKNESS)
    outer_boundary = openmc.Sphere(r=OUTER_BOUNDARY_R, boundary_type="vacuum")

    void_inner = openmc.Cell(
        name="30: void outside shield, inside detector shell",
        region=+shell_surfaces[-1] & -detector_inner,
    )
    detector_shell = openmc.Cell(
        name="F2-equivalent tally shell at r=160 cm",
        region=+detector_inner & -detector_outer,
    )
    void_outer = openmc.Cell(
        name="40: void past detector",
        region=+detector_outer & -outer_boundary,
    )

    model.geometry = openmc.Geometry(
        [cavity, *shield_cells, void_inner, detector_shell, void_outer]
    )

    # --- Source: point isotropic 7-MeV photon source -----------------
    source = openmc.IndependentSource()
    source.particle = "photon"
    source.space = openmc.stats.Point((0.0, 0.0, 0.0))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([SOURCE_ENERGY_EV], [1.0])
    model.settings.source = source

    # --- Settings ------------------------------------------------------
    model.settings.run_mode = "fixed source"
    model.settings.photon_transport = True
    model.settings.electron_treatment = "led"  # no bremsstrahlung, cf. phys:p
    model.settings.batches = 20
    model.settings.particles = 50_000

    # --- Tallies ---------------------------------------------------
    detector_filter = openmc.CellFilter(detector_shell)
    dose_filter = openmc.EnergyFunctionFilter(DE2_MEV * 1.0e6, DF2_SV_CM2)

    flux_tally = openmc.Tally(name="surface flux at r=160 cm")
    flux_tally.filters = [detector_filter]
    flux_tally.scores = ["flux"]

    dose_tally = openmc.Tally(name="ambient dose equivalent H*(10) at r=160 cm")
    dose_tally.filters = [detector_filter, dose_filter]
    dose_tally.scores = ["flux"]

    model.tallies = openmc.Tallies([flux_tally, dose_tally])

    model._detector_shell_volume = (
        4.0 / 3.0 * np.pi * (detector_outer.r**3 - detector_inner.r**3)
    )
    return model


def run_and_report(include_shield: bool = True, **run_kwargs):
    model = build_model(include_shield=include_shield)
    volume = model._detector_shell_volume

    sp_path = model.run(**run_kwargs)
    with openmc.StatePoint(sp_path) as sp:
        dose = sp.get_tally(name="ambient dose equivalent H*(10) at r=160 cm")
        mean = dose.mean.flatten()[0] / volume
        std_dev = dose.std_dev.flatten()[0] / volume

    label = "WITH iron shield" if include_shield else "bare source (no shield)"
    print(f"\n=== {label} ===")
    print(f"Ambient dose equivalent H*(10) at r=160 cm: "
          f"{mean:.4e} +/- {std_dev:.4e} Sv/photon")
    if not include_shield:
        print("Book reference (no shield): 6.013e-17 Sv/photon")
    return mean, std_dev


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-shield", action="store_true",
        help="Run the bare-source validation case instead of the shielded case",
    )
    args = parser.parse_args()
    run_and_report(include_shield=not args.no_shield)
