Simulation Control
==================

This defines the settings of the simulation period and timestep.

**Example:**

.. code-block:: yaml

    Simulation Control:
      Start Date: 2015-Jan-1
			End Date: 2015-Dec-31
			Timestep: 60 # [min]

Start Date
----------

Specifies the start date of the simulation. Simulation begins at 12:00am of this day.

=============   =========================================
**Required:**   Yes
**Type:**       Date String (eg: YYYY-Mon-DD, YYYY/MM/DD)
=============   =========================================

End Date
--------

Specifies the end date of the simulation. Simulation ends one timestep previous to 12:00am of the following day.

=============   =========================================
**Required:**   Yes
**Type:**       Date String (eg: YYYY-Mon-DD, YYYY/MM/DD)
=============   =========================================

Timestep
--------

Timestep duration in minutes used in calculations.

=============   =======
**Required:**   Yes
**Type:**       Integer
**Units:**      Minutes
=============   =======
