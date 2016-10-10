Initialization
==============

Defines how the inital temperatures are generated for the model.

**Example:**

.. code-block:: yaml

  Initialization:
    Initialization Method: STEADY-STATE  # KUSUDA | CONSTANT | STEADY-STATE
    Accelerated Initialization Timestep: 168  # hours
    Number of Accelearted Initialization Timesteps: 12
    Number of Warmup Days in Initialization: 365 # days

Initialization Method
---------------------

The initialization method determines how the initial temperatures in the domain are set. Options are:

- ``CONSTANT``, spatially-constant initial temperature,
- ``KUSUDA``, a one-dimensional analytical solution developed by that provides temperature variation as a function of depth driven by an annual harmonic temperature fluctuation. There is no temperature variation in horizontal dimensions,
- ``STEADY-STATE``, a steady-state solution scheme initializes the temperatures with the first timestep’s boundary conditions. This provides an initial condition temperature variation in all dimensions.

=============   =============================================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``CONSTANT``, ``KUSUDA``, or ``STEADY-STATE``
**Default:**    ``STEADY-STATE``
=============   =============================================

Initial Temperature
-------------------

When `Initialization Method`_ is ``CONSTANT`` this specifies the temperature to use.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

Accelerated Initialization Timestep
-----------------------------------

An accelerated initialization begins with the user-defined `Initialization Method`_ and calculates new domain temperatuers prior to the beginning of the simulation using long timesteps (on the order of days, weeks, or months). These timesteps are calculated using a fully implicit, unconditionally stable numerical scheme. This allows the simulation to build a history of temperatures without requiring a signficant amount of additional calculations. The defualt, one week, was found to give very accurate initial temperatures.

=============   =======
**Required:**   No
**Type:**       Integer
**Units:**      hours
**Default:**    168
=============   =======

Number of Accelerated Initialization Timesteps
----------------------------------------------

This specifies the number of timesteps (of the size specified by `Accelerated Initialization Timestep`_) to calculate prior to the beginning of the simulation.

=============   =======
**Required:**   No
**Type:**       Integer
**Default:**    12
=============   =======

Number of Warmup Days in Initialization
---------------------------------------

Additional days of initialization can be calculated using the :ref:`timestep` and `Numerical Scheme`_ defined by the user. This input specifies the number of days the domain is simulated under these conditions after the accelerated initialization timesteps, but prior to the :ref:`start_date` specified in the :ref:`simulation_control`.

=============   =======
**Required:**   No
**Type:**       Integer
**Units:**      days
**Default:**    365
=============   =======
