Data Sources
============

pynucastro makes use of data sources from several different nuclear
physics compilations.  Here we explain how these are used in
pynucastro.

ReacLib
-------

The `JINA ReacLib database <https://reaclib.jinaweb.org/>`_ provides
the core thermonuclear rates.  These are managed in pynucastro through
the :func:`ReacLibRate <pynucastro.rates.rate.ReacLibRate>` class.

A ReacLib rate is composed of several sets (managed as
:func:`ReacLibRate <pynucastro.rates.rate.SingleSet>` objects, each
taking the form:

.. math::

   \lambda = \exp{\left (a_0 + \sum_{i=1}^5  a_i T_9^{(2i-5)/3}  + a_6 \log T_9\right )}

where :math:`\lambda` is of the form:

.. math::

   \lambda = \left \{ \begin{array}{cc} \tau^{-1}  & \mbox{unary reaction} \\
                                     N_A \langle \sigma v\rangle & \mbox{binary reaction} \\
                                     N_A^2 \langle \sigma v\rangle & \mbox{trinary reaction}
                   \end{array} \right .


The rates are divided into chapters, based on how many nuclei are on each side of the reaction:

========  ====================================================
chapter    reaction format
--------  ----------------------------------------------------
1         :math:`e_1 \rightarrow e_2`
2         :math:`e_1 \rightarrow e_2 + e_3`
3         :math:`e_1 \rightarrow e_2 + e_3 + e_4`
4         :math:`e_1 + e_2 \rightarrow e_3`
5         :math:`e_1 + e_2 \rightarrow e_3 + e_4`
6         :math:`e_1 + e_2 \rightarrow e_3 + e_4 + e_5`
7         :math:`e_1 + e_2 \rightarrow e_3 + e_4 + e_5 + e_6`
8         :math:`e_1 + e_2 + e_3 \rightarrow e_4`
9         :math:`e_1 + e_2 + e_3 \rightarrow e_4 + e_5`
10        :math:`e_1 + e_2 + e_3 + e_4 \rightarrow e_5 + e_6`
11        :math:`e_1 \rightarrow e_2 + e_3 + e_4 + e_5`
========  ====================================================


Sources

Electron capture

Reverse rates


