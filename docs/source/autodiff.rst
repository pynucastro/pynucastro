Automatic differentiation
=========================

The :py:mod:`autodiff <pynucastro.autodiff>` module provides basic
forward-mode automatic differentiation routines for python functions.

.. note::

   For Numba compatibility, these are implemented as a python
   ``namedtuple`` and not a class.  Numba's ``jitclass`` does
   not associate reflected operators (i.e. scalar + a ``Dual`` object)
   with a class.  The ``namedtuple`` approach was recommended by
   the Numba developers (see `Numba issue #10921 <https://github.com/numba/numba/issues/10921>`_).

.. note::

   For :py:obj:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>`, we use the C++ `autodiff library <https://github.com/autodiff/autodiff>`_.


Single derivative
-----------------

:py:obj:`Dual <pynucastro.autodiff.dual.Dual>`
provides forward-mode automatic differentiation with one
derivative per value.

Here's an example.  We seed the independent variable with ``Dual(value,
1.0)`` (the default derivative is one):

.. code:: python

   from numba import njit
   from pynucastro.autodiff.dual import Dual

   @njit
   def polynomial(x):
       return x * x + 3.0 * x - 1.0

   result = polynomial(Dual(1.3))

This gives:

.. code:: python

   Dual(value=4.590000000000001, deriv=5.6)

``Dual`` is an immutable namedtuple with fields ``value`` and ``deriv``.
Numba operator overloads support dual-dual arithmetic and real scalars on
either side of addition, subtraction, multiplication, and division.
Powers support a dual base and a scalar or dual exponent; a dual exponent
requires a positive base.

.. important::

   For pure python, the standard operator functions, e.g., ``__add__()``
   are invoked for arithmetic.  When compiled with Numba, separate
   functions decorated with ``@overload`` provide the implementations.
   Both use the same derivative formulas.  Numba represents the ``namedtuple``
   as a native structure rather than a heap-allocated ``jitclass``.

.. tip::

   Constants can simply be written as scalars or explicitly as a
   ``Dual`` with a seed of zero, i.e., ``Dual(value, 0.0)``.


Math functions
^^^^^^^^^^^^^^

The module's version of math functions, ``exp``, ``log``, ``log1p``,
``sqrt``, ``cbrt``, ``sin``, ``cos``, ``arctan``, should be used when
an argument may be a dual.  These fall-back to their Numpy equivalents
for scalar input.  The derivative of ``abs`` is defined as zero at zero.

Comparisons (e.g. ``<``, ``>``) use only the value.

Compiled two-argument ``min`` and ``max`` return a dual with a
floating-point value, even when a scalar bound wins; that bound has
zero derivative components with the same seed kind as the dual
operand.  Ties select the first operand. Python's built-in ``min`` and
``max`` return the selected operand itself, which may be a scalar. In
other compiled branches, keep the return types compatible explicitly.

Two independent derivatives
===========================

:py:obj:`DerivativePair
<pynucastro.autodiff.derivative_pair.DerivativePair>` differentiates
with respect to two variables in one evaluation. Both components are
scalar fields in a namedtuple, so this does not require a NumPy
derivative array::

    from pynucastro.autodiff.derivative_pair import DerivativePair

    x = Dual(2.0, DerivativePair(1.0, 0.0))
    y = Dual(3.0, DerivativePair(0.0, 1.0))
    result = x * y
    # result.value = 6; result.deriv.d1 = 3; result.deriv.d2 = 2

All dual operands in an expression must use the same kind of seed.
Constants may be ordinary real scalars. ``constant_like(x, value)`` creates
a constant with the same scalar or dual representation as ``x``, including
zero derivative components; this is useful for compiled clipping branches.
