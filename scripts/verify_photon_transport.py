#!/usr/bin/env python3
"""Smoke test: run a fixed-source photon transport calculation through a
lead sphere and print the results.

Requires the OPENMC_CROSS_SECTIONS environment variable to point at a
cross_sections.xml built by setup_nuclear_data.sh.
"""
import tempfile
from pathlib import Path

import openmc


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        prev_dir = Path.cwd()
        try:
            import os

            os.chdir(tmpdir)

            lead = openmc.Material(name="lead")
            lead.add_element("Pb", 1.0)
            lead.set_density("g/cm3", 11.35)
            openmc.Materials([lead]).export_to_xml()

            sphere = openmc.Sphere(r=10.0, boundary_type="vacuum")
            cell = openmc.Cell(fill=lead, region=-sphere)
            openmc.Geometry([cell]).export_to_xml()

            settings = openmc.Settings()
            settings.run_mode = "fixed source"
            settings.photon_transport = True
            settings.particles = 1000
            settings.batches = 5
            settings.source = openmc.IndependentSource(
                particle="photon",
                space=openmc.stats.Point((0, 0, 0)),
                energy=openmc.stats.Discrete([1.0e6], [1.0]),
            )
            settings.export_to_xml()

            tally = openmc.Tally(name="flux")
            tally.scores = ["flux", "heating"]
            openmc.Tallies([tally]).export_to_xml()

            openmc.run()

            with openmc.StatePoint("statepoint.5.h5") as sp:
                t = sp.get_tally(name="flux")
                print(t.get_pandas_dataframe())
        finally:
            os.chdir(prev_dir)


if __name__ == "__main__":
    main()
