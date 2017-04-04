Boundaries
==========

Definitions related to the boundary conditions of the ground domain.

**Example:**

.. code-block:: yaml

  Boundaries:
    Far-Field Width: 40.0  # [m]
    Deep-Ground Depth: 40.0  # [m]
    Deep-Ground Boundary Condition: ZERO-FLUX  # AUTO | CONSTANT-TEMP | ZERO-FLUX
    Indoor Air Temperature: 295.0  # [K]

Indoor Air Temperature Method
-----------------------------

Allows the user to choose between having a constant indoor temperature for the duration of the simulaiton or to reference temperatures from a file.

=============   ========================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``FILE`` or ``CONSTANT``
**Default:**    ``CONSTANT``
=============   ========================

Indoor Air Temperature File
---------------------------

If `Indoor Air Temperature Method`_ is ``FILE`` the indoor dry-bulb temperature (in Kelvin) will be set using hourly values defined in a comma separted value (CSV) file.

**Example:**

.. code-block:: yaml

  Indoor Air Temperature File:
    Name: ../path/to/file.csv
    Index: [1,1]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============


Name
^^^^

Path (relative or absolute) file.

=============   =========
**Required:**   Yes
**Type:**       File Path
=============   =========

Index
^^^^^

A list of two values corresponding to the row and column where the hourly data begins in the file. A value of ``[0, 0]`` starts at the first row and first column. A value of ``[0,1]`` starts at the first row and second column.

=============   ====================
**Required:**   Yes
**Type:**       List [2] of integers
=============   ====================

.. _in_temp:

Indoor Air Temperature
----------------------

If `Indoor Air Temperature Method`_ is ``CONSTANT`` the indoor dry-bulb temperature will be set using this value. If `Indoor Air Temperature Method`_ is ``FILE``, then this is not required.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

Outdoor Air Temperature Method
------------------------------

Allows the user to choose between having a constant outdoor temperature for the duration of the simulaiton or to reference temperatures from the weather file.

=============   ================================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``WEATHER-FILE`` or ``CONSTANT``
**Default:**    ``WEATHER-FILE``
=============   ================================

.. _out_temp:

Outdoor Air Temperature
-----------------------

If `Outdoor Air Temperature Method`_ is ``CONSTANT`` the outdoor dry-bulb temperature will be set using this value. If `Outdoor Air Temperature Method`_ is ``WEATHER-FILE``, then this is not required.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

Local Boundary Layer Thickness
------------------------------

Local boundary layer thickness used for calculating local wind speeds from weather file wind speeds.

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      m
**Default:**    370
=============   =======

Local Terrain Exponent
----------------------

Local terrain exponent used for calculating local wind speeds from weather file wind speeds.

=============   =============
**Required:**   No
**Type:**       Numeric
**Units:**      dimensionless
**Default:**    0.22
=============   =============

Far-Field Width
---------------

Distance from the interior wall surface to the edge of the domain.

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      m
**Default:**    40
=============   =======

Deep-Ground Depth
-----------------

Distance from the grade level to the bottom of the domain.

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      m
**Default:**    40
=============   =======

Deep-Ground Boundary Condition
------------------------------

Specifies the type of boundary condition to apply at the deep-ground boundary. Options are:

- ``ZERO-FLUX``, which applies a zero heat flux boundary,
- ``AUTO``, which applies a constant temperature equal to the average outdoor dry-bulb temperature from the weather file, and
- ``CONSTANT-TEMP``, which applies a user-specified constant temperature (see `Deep-Ground Temperature`_).

=============   =============================================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``ZERO-FLUX``, ``AUTO``, or ``CONSTANT-TEMP``
**Default:**    ``ZERO-FLUX``
=============   =============================================

Deep-Ground Temperature
-----------------------

If `Deep-Ground Boundary Condition`_ is ``CONSTANT-TEMP``, then this value specifies the temperature applied to the deep-ground boundary.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

Convection Calculation Method
-----------------------------

Specifies how convection coefficients are calculated. Options are:

- ``AUTO``, which calculates dynamic convection coefficients based on temperature difference, wind speed, and wind direction.
- ``CONSTANT``, which applies a user-specified convection coefficients to interior and exterior surfaces (see `Interior Convection Coefficient`_ and `Exterior Convection Coefficient`_). This is used primariliy for IEA BESTEST calculations.

=============   ========================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``AUTO`` or ``CONSTANT``
**Default:**    ``AUTO``
=============   ========================

Interior Convection Coefficient
-------------------------------

If `Convection Calculation Method`_ is ``CONSTANT``, then this value specifies the convection coefficient applied to interior surface boundaries (slab floor, interior foundation wall, and interior insulation).

=============   ===============
**Required:**   Depends
**Type:**       Numeric
**Units:**      W/m\ :sup:`2`-K
=============   ===============

Exterior Convection Coefficient
-------------------------------

If `Convection Calculation Method`_ is ``CONSTANT``, then this value specifies the convection coefficient applied to exterior surface boundaries (grade, exterior foundation wall, and exterior insulation).

=============   ===============
**Required:**   Depends
**Type:**       Numeric
**Units:**      W/m\ :sup:`2`-K
=============   ===============

Wall Top Boundary Condition
---------------------------

Specifies how the boundary condition along the wall top is calculated. Options are:

- ``ZERO-FLUX``, which applies a zero heat flux boundary condition along the wall top. This implies that heat flux above the wall top is one dimensional and does not flow through the wall top boundary.
- ``LINEAR-DT``, which applies a linear change in temperature across the wall top (see `Wall Top Temperature Difference`_). This is used only to represent the constraints of the IEA BESTEST analytical solution in case GC10a.

=============   ==============================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``ZERO-FLUX`` or ``LINEAR-DT``
**Default:**    ``ZERO-FLUX``
=============   ==============================

Wall Top Interior Temperature
-----------------------------

If `Wall Top Boundary Condition`_ is ``LINEAR-DT``, then this value specifies the interior temperature at the wall top. This is used only to represent the constraints of the IEA BESTEST analytical solution in case GC10a.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

Wall Top Exterior Temperature
-----------------------------

If `Wall Top Boundary Condition`_ is ``LINEAR-DT``, then this value specifies the exterior temperature at the wall top. This is used only to represent the constraints of the IEA BESTEST analytical solution in case GC10a.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======
