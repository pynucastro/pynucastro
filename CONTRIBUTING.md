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
done via pull requests.* A pull request should be generated from your
fork of pynucastro and target the `main` branch.

You should run the `py.test` unit tests on your changes, if possible,
before issuing the PR. To run the unit tests, in the top `pynucastro/`
directory, run:

```
py.test-3 -v .
```

Once you have run the unit tests and submitted the PR, one of the
pynucastro developers will review the PR and if needed, suggest
modifications prior to merging the PR.

If there are a number of small commits making up the PR, we may wish
to squash commits upon merge to have a clean history.  *Please ensure
that your PR title and first post are descriptive, since these will be
used for a squashed commit message.*
