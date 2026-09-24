#!/usr/bin/env bash
# Downloads ENDF/B-VIII.0 (NNDC-processed) nuclear data and builds a
# cross_sections.xml suitable for OpenMC photon transport.
#
# OpenMC's photon transport requires BOTH:
#   - incident photon (photoatomic) data, one HDF5 file per element
#   - incident neutron data, one HDF5 file per isotope (OpenMC uses this
#     to resolve atomic weight ratios for every nuclide referenced by a
#     material, even in photon-only runs)
#
# Data source: https://github.com/openmc-data-storage/ENDF-B-VIII.0-NNDC
# (a pre-converted HDF5 mirror of the official NNDC ENDF/B-VIII.0 release,
# used because the primary NNDC/openmc.org hosts may not be reachable from
# all network environments).
#
# Usage: ./setup_nuclear_data.sh [destination_dir]
# Default destination: $HOME/nuclear_data

set -euo pipefail

DEST="${1:-$HOME/nuclear_data}"
CLONE_DIR="$(mktemp -d)"
REPO_URL="https://github.com/openmc-data-storage/endf-b-viii.0-nndc"

trap 'rm -rf "$CLONE_DIR"' EXIT

echo "Cloning nuclear data repository (this is several GB, may take a while)..."
GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 "$REPO_URL" "$CLONE_DIR"

mkdir -p "$DEST/neutron" "$DEST/photon"

echo "Copying photon (photoatomic) data..."
cp "$CLONE_DIR"/h5_files/photon/*.h5 "$DEST/photon/"

echo "Copying neutron data (skipping any unresolved Git LFS pointer files)..."
find "$CLONE_DIR"/h5_files/neutron -name '*.h5' -size +1000c -exec cp {} "$DEST/neutron/" \;

missing=$(find "$CLONE_DIR"/h5_files/neutron -name '*.h5' -size -1000c | wc -l)
if [ "$missing" -gt 0 ]; then
    echo "NOTE: $missing neutron file(s) were Git LFS pointers (not fetched by the" >&2
    echo "anonymous clone) and were skipped. Materials using those isotopes" >&2
    echo "(e.g. U235, U238) will not be usable until real data is obtained." >&2
fi

echo "Building cross_sections.xml..."
python3 - "$DEST" <<'PYEOF'
import sys
from pathlib import Path
import openmc.data

dest = Path(sys.argv[1])
lib = openmc.data.DataLibrary()
for h5file in sorted((dest / 'neutron').glob('*.h5')):
    lib.register_file(h5file)
for h5file in sorted((dest / 'photon').glob('*.h5')):
    lib.register_file(h5file)
lib.export_to_xml(dest / 'cross_sections.xml')
print(f"Wrote {dest / 'cross_sections.xml'}")
PYEOF

echo
echo "Done. Set the following environment variable before running OpenMC:"
echo "  export OPENMC_CROSS_SECTIONS=$DEST/cross_sections.xml"
