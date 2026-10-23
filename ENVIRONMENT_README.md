# Software environments

Two Python environments are provided because the original pipeline and the revision
analyses were run with different software.

| File | Scripts | Status |
|---|---|---|
| `environment_original.yml` | `Source Code with Raw Data/*.py` | Versions taken from script headers |
| `environment_updated.yml` | `Updated Analysis/**/*.py` | Only PyDESeq2 0.5.4 documented; remaining versions to be filled in |

Both were install- and import-tested (see header comments). The original environment needs `setuptools<81`, and the updated one needs Python 3.11 or later because of PyDESeq2 0.5.4.

Create an environment with `conda env create -f environment_original.yml`.

R scripts record their R and package versions in their header comments
(R 4.3.3 to 4.4.2; see each script).
