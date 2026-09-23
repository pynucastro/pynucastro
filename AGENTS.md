# Coding-agent guide for pynucastro

## Project purpose and layout

`pynucastro` is a Python library for nuclear astrophysics: nuclear data,
reaction-rate libraries, reaction-network construction and analysis, and code
generation for Python, C++, Fortran, and AMReX-Astro backends.  Numerical and
scientific correctness is part of the public API: do not make a seemingly
cosmetic change to rate evaluation, units, reaction identities, or generated
code without targeted regression tests.

Key locations:

- `pynucastro/`: library source.  The main domains are `rates/`, `nucdata/`,
  `networks/`, `screening/`, `eos/`, `neutrino_cooling/`, and `reduction/`.
- `pynucastro/<domain>/tests/`: unit and regression tests colocated with each
  domain.  Network-generation and backend-comparison tests are under
  `pynucastro/networks/tests/`.
- `pynucastro/data/` and `pynucastro/nucdata/`: packaged scientific data.
  Preserve source/provenance, file formats, and units when changing them.
- `pynucastro/templates/`: templates used to generate exported networks.  A
  template change normally requires regenerating and checking reference
  networks.
- `docs/source/`: Sphinx documentation (reStructuredText and MyST notebooks).
  `docs/source/index.rst` is the documentation table of contents and
  `docs/source/conf.py` is the build configuration.
- `.github/workflows/`: the CI commands and platform/backend coverage.

## Public API map

`pynucastro/__init__.py` is the authoritative high-level map of the public
classes.  Consult its module docstring before changing a public interface.
The following summary is intentionally compact; inspect the defining module
and its tests for implementation details.

- Nuclear data: `Nucleus` represents a nuclide and is the primary interface to
  mass, spin, and partition-function data; `Composition` stores abundances for
  a set of nuclei.
- Rates and libraries: `Rate` is the base reaction-rate abstraction.
  `ReacLibRate`, `TabularWeakRate`, `TemperatureTabularRate`, and `StarLibRate`
  provide the principal source/representation-specific evaluations.  `Library`
  stores and filters rates, with `RateFilter` defining queries and specialized
  libraries loading the supported datasets.
- Rate transformations: `DerivedRate` is a detailed-balance reverse rate;
  `ModifiedRate` changes an underlying rate's representation/stoichiometry;
  `BranchedRate` models endpoint branching; and `ApproximateRate` computes an
  effective rate from hidden rates.  These are scientifically meaningful
  wrappers, not interchangeable conveniences.
- Network core: `RateCollection` organizes linked nuclei and rates and
  evaluates or visualizes a network.  `Explorer` provides interactive notebook
  exploration, while `NSENetwork` extends the collection with nuclear
  statistical-equilibrium solving.
- Exporting networks: `PythonNetwork`, `SimpleCxxNetwork`, and
  `AmrexAstroCxxNetwork` generate ODE right-hand-side implementations for
  their respective targets.  `BaseCxxNetwork` is their C++ base class;
  `FortranNetwork` provides a Fortran interface around the simple-C++ path.
- Other physics: `StellarEOS` and `FermiIntegral` support the equation of
  state; `NeutrinoCooling` models thermal neutrino losses; screening helpers
  modify charged-particle rates; and `reduction` provides DRGEP and sensitivity
  analysis tools.

## Working conventions

- Work from the repository root unless a command below states otherwise.
- Before editing, inspect `git status --short`.  This repository may contain
  user-created notebooks, generated networks, compiled executables, and other
  artifacts; never delete, revert, or fold unrelated changes into a patch.
- Keep changes scoped.  Add or update a test whenever behavior changes; add
  user-facing documentation for public API, scientific behavior, or workflow
  changes.
- Preserve established numerical behavior unless the task deliberately changes
  it.  State the expected physical/numerical effect, units, and tolerance in
  tests or documentation where that is useful.
- Use the existing `Rate` identity conventions.  `Rate.id` and especially
  programming-safe `Rate.fname` must remain unique in libraries and generated
  networks.  Copy a rate before mutating a rate obtained from a `Library`.
- Do not hand-edit generated network reference outputs as a substitute for
  fixing generation code.  Regenerate them only when the intended output has
  changed, then review the diff.

## Environment and tests

Install the development and documentation dependencies, then install the
package in editable mode:

```bash
python -m pip install -r requirements.txt -r requirements-docs.txt
python -m pip install --editable .
```

Use the narrowest relevant test first:

```bash
pytest -v pynucastro/rates/tests
pytest -v pynucastro/networks/tests
pytest -v                         # unit tests
pytest -v --nbval                 # unit tests plus notebooks
```

`conftest.py` configures NumPy CPU features before NumPy is imported so that
network-output regression tests are reproducible.  Run tests through `pytest`
rather than bypassing that setup.

Changes to network writers or C++/Fortran/AMReX templates may require:

```bash
pytest -v -s --update-networks -k write_network
```

This command intentionally updates reference outputs.  Run it only for an
expected generation change and inspect every changed artifact.  Backend tests
may also require local compilers or AMReX/Microphysics dependencies; do not
mask a missing optional toolchain by weakening a test.

## Style and static checks

Configuration lives in `pyproject.toml`.  CI runs ruff, pylint, flake8,
pydocstyle, isort, codespell, pytest, and backend-specific workflow tests.
For changed Python code, run the applicable checks when available:

```bash
ruff check pynucastro
isort --check --diff pynucastro
pylint pynucastro pynucastro/**/tests
flake8 pynucastro
pydocstyle pynucastro
```

Follow the existing local style rather than reformatting unrelated files.
Python line length is 132; public functions should have NumPy-style docstrings
(the project has a small, documented set of pydocstyle exceptions).

## Documentation

For documentation changes, build from `docs/`:

```bash
cd docs
SKIP_EXECUTE=TRUE make html
```

Without `SKIP_EXECUTE=TRUE`, the configured MyST notebooks execute during the
build.  The docs Makefile also downloads the current Zenodo citation and
generates `docs/source/changelog.md` from `CHANGES.md`; these are generated
files, so review their diff and do not edit them manually.  CI uses the
stricter build:

```bash
make SPHINXOPTS='-v -W --keep-going -n' html
make SPHINXOPTS=-v linkcheck
```

If an intentionally changed notebook output needs recording, use
`docs/regen_notebook.sh <notebook paths>` and review the notebook diff.

## Completion checklist

- Run focused tests and relevant static checks; report any checks not run and
  why.
- For answer-changing changes, call out the scientific effect and regression
  coverage.
- For templates or generated outputs, confirm the regenerated files are
  intentional.
- Leave unrelated working-tree changes untouched.
