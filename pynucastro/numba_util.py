"""Drop-in for numba that provides no-op stubs if numba isn't installed."""

try:
    # the linters really don't like this
    # pylint: disable=unused-import, wildcard-import, unused-wildcard-import
    # flake8: noqa: F401, F403
    from numba import jit, njit
    from numba.experimental import jitclass
    from numba.types import *
except ImportError:
    import functools

    def _noop_wrapper(obj):
        # this copies the function name and docstring to the wrapper function
        @functools.wraps(obj)
        def wrapper(*args, **kwargs):
            return obj(*args, **kwargs)
        return wrapper

    def jit(*args, **kwargs):
        """Apply just-in-time compilation."""
        if len(args) == 1 and not kwargs and callable(args[0]):
            # used as a normal decorator with no arguments (`@jit`)
            return _noop_wrapper(args[0])
        # used as a decorator factory (`@jit(...)`)
        return _noop_wrapper

    njit = jit

    def jitclass(*args, **kwargs):
        """Use just-in-time compilation for a class."""
        if len(args) == 1 and not kwargs and isinstance(args[0], type):
            # used as a normal decorator
            return _noop_wrapper(args[0])
        # used as a decorator factory
        return _noop_wrapper

    # placeholders for numba types, needed for jitclass specs
    class _FakeNumbaType:
        def __getitem__(self, _key):
            return self

    # assume everything not defined above is a type
    def __getattr__(_name):
        return _FakeNumbaType()
