.. _simulation_control:

Simulation Control
==================

This defines the settings of the simulation period and timestep.

**Example:**

.. code-block:: yaml

    Simulation Control:
      Start Date: 2015-Jan-1
      End Date: 2015-Dec-31
      Timestep: 60 # [min]

.. _start_date:

Start Date
----------

Specifies the start date of the simulation. Simulation begins at 12:00am of this day. This is specified as a date string (e.g., YYYY-Mon-DD, YYYY/MM/DD).

=============   ===========
**Required:**   Yes
**Type:**       Date string
=============   ===========

End Date
--------

Specifies the end date of the simulation. Simulation ends one timestep previous to 12:00am of the following day. This is specified as a date string (e.g., YYYY-Mon-DD, YYYY/MM/DD).

=============   ===========
**Required:**   Yes
**Type:**       Date string
=============   ===========

.. _timestep:

Timestep
--------

Timestep duration in minutes used in calculations.

=============   =======
**Required:**   Yes
**Type:**       Integer
**Units:**      Minutes
=============   =======
