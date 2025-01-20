# Contributing

Contributions are welcomed from anyone, including posting issues or
submitting pull requests to the [pynucastro github](https://github.com/pynucastro/pynucastro).

## Issues

Creating an issue on github is a good way to request new features,
file a bug report, or notify us of any difficulties that arise using
pynucastro.

To request support using pynucastro, please create an issue on the
pynucastro github and the developers will be happy to assist you.

If you are reporting a bug, please indicate any information necessary to
reproduce the bug including your version of python.

## Pull Requests

*Any contributions that have the potential to change answers should be
done via pull requests.* A pull request should be generated from a
feature branch in your fork of pynucastro and target the `main` branch.

You should run the test suite on your changes, if possible, before
issuing the PR. See the section below for detailed instructions.

Once you have run the tests and submitted the PR, one of the
pynucastro developers will review the PR and if needed, suggest
modifications prior to merging the PR.

If there are a number of small commits making up the PR, we may wish
to squash commits upon merge to have a clean history.  *Please ensure
that your PR title and first post are descriptive, since these will be
used for a squashed commit message.*

## Documentation

The documentation is in `docs/source` and uses MyST-flavored markdown,
managed by Sphinx.  To build the docs, in the `docs/` directory, do:

```
make html
```

by default, this will execute all of the Jupyter notebooks.  To skip this,
build as:

```
SKIP_EXECUTE=TRUE make html
```

## Testing

We use pytest to do unit and regression tests. All commands should be
run from the repository root.

* To run all the tests:

  ```
  pytest -v --nbval
  ```

* To run just the unit tests (this is faster, but may not catch all bugs):

  ```
  pytest -v
  ```

* To run just the jupyter notebook regression tests:

  ```
  pytest -v --nbval -p no:python
  ```

  This executes each of the notebooks under `examples/` and `docs/source/`
  and compares against the stored outputs from each cell.

* To check code coverage:

  ```
  pytest --cov=pynucastro --nbval
  ```

  The results can be inspected in more detail by running `coverage html`
  and opening the generated `htmlcov/index.html` in a web browser.

  Note that numba-accelerated routines (most notably the screening
  functions) do not support coverage reporting, and will always show up
  as missed.

* To regenerate the reference network output files used in the unit tests:

  ```
  pytest -v -s --update-networks -k write_network
  ```

  This will be needed if any of the network generation code or C++
  network templates are changed.

* To re-run notebooks whose outputs have changed:

  ```
  docs/regen_notebook.sh <paths to notebooks>
  ```
