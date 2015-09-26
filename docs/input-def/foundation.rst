Foundation
==========

The description of the foundation design is provided within the two-dimensional context. This profile is applied along the entire perimeter of the foundation.

.. figure:: ../images/context.png

   Two-dimensional context for ``Foundation`` object definition

The foundation insulation and structural components are defined by (up to) six sub-objects. These objects allow the user to flexibly describe any foundation construction.

.. figure:: ../images/components.png

   Insulation and structural design components

**Example:**

.. code-block:: yaml

  Foundation:
    Foundation Depth: 0.0 # [m]
    Polygon:
      - [0, 0] # [m, m]
      - [0, 20] # [m, m]
      - [20, 20] # [m, m]
      - [20, 0] # [m, m]
    Soil: Typical Soil # Material reference
    Slab:
      Layers:
        -
          Material: Concrete # Material reference
          Thickness: 0.2032 # [m]
    Wall:
      Layers:
        -
          Material: Concrete # Material reference
          Thickness: 0.3048 # [m]
      Height Above Grade: 0.3048  # [m]
      Height: 0.508 # [m]
    Interior Horizontal Insulation:
      Depth: 0.2032 # [m]
      Width: 0.4064 # [m]
      Material: XPS # Material reference
      Thickness: 0.0508
    Interior Vertical Insulation:
      Depth: 0.2032 # [m]
      Material: XPS # Material reference
      Thickness: 0.0508 # [m]
    Indoor Air Temperature: 295.372 # [K]


Foundation Depth
----------------

``Foundation Depth`` defines the distance from the wall top to the floor. This value is used to characterize the type of foundation (slab, crawlspace, or basement). For example, a value of zero would represent a sla and a value near 2 meters would represent a basement.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Polygon
-------

The foundation shape is defined by the description of a single polygon. The perimeter of this polygon defines the location of the interior surface of the foundation wall. The positioning of the foundation insulation and structural components are translated into three-dimensional space internally.

The polygon is specified by a list of x-y Cartesian vertices tracing the foundation perimeter in a clockwise fashion. When simulating in three-dimensions, the polygon must be rectilinear (comprised only of right angles).

.. figure:: ../images/polygon.png

   Plan view illustrating foundation shape vertex definition and far-field boundary.

**Example:**

.. code-block:: yaml

    Polygon:
      - [0, 20]
      - [15, 20]
      - [15, 30]
      - [30, 30]
      - [30, 17]
      - [22, 17]
      - [22, 0]
      - [12, 0]
      - [12, 10]
      - [0, 10]

=============   =================================
**Required:**   Yes
**Type:**       List [N] of lists [2] of numerics
**Units:**      m
=============   =================================


Soil
----

Represents the soil surrounding the building foundation.

=============   ==================
**Required:**   Yes
**Type:**       Material reference
=============   ==================


Slab
----

This defines the costruction of the floor slab in the foundation. This is not required. If there is no slab defined for a given foundation, then the floor is exposed soil.

**Example:**

.. code-block:: yaml

  Slab:
    Layers:
      -
        Material: XPS # Material reference
        Thickness: 0.0508 # [m]
      -
        Material: Concrete # Material reference
        Thickness: 0.2032 # [m]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============

Layers
^^^^^^

Layers are specified as a list of material references, and thickness values from the outtermost layer to the innermost layer (at the floor surface). A layer of insulation can be added to model whole-slab insulation.

=============   ===========================================
**Required:**   Yes
**Type:**       List of layers (a material and a thickness)
=============   ===========================================

Material
""""""""

Material composing the layer.

=============   ==================
**Required:**   Yes
**Type:**       Material reference
=============   ==================

Thickness
"""""""""

Thickness of the layer.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======


Emissivity
^^^^^^^^^^

Interior emissivity of the slab used for interior long-wave radiation calculations.

=============   =============
**Required:**   No
**Type:**       Numeric
**Units:**      dimensionless
**Default:**    0.8
=============   =============

Wall
----

This defines the costruction of the foundation wall. This is not required. If there is no wall defined for a given foundation, then the wall is exposed soil.

**Example:**

.. code-block:: yaml

  Wall:
    Height: 2.95 # [m]
    Height Above Grade: 0.3048  # [m]
    Layers:
      -
        Material: XPS # Material reference
        Thickness: 0.0508 # [m]
      -
        Material: Concrete # Material reference
        Thickness: 0.2032 # [m]
      -
        Material: XPS # Material reference
        Thickness: 0.0508 # [m]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============

Height
^^^^^^

The height of the wall describes the distance from the wall top to the bottom of the foundation footer (the footer is not modeled separately). This value should generally be greater than that of the `Foundation Depth`_ combined with the total thickness of the slab.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======


Height Above Grade
^^^^^^^^^^^^^^^^^^

The height of the wall top relative to the grade (z = 0).

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======


Layers
^^^^^^

Layers are specified as a list of material references, and thickness values from the outtermost layer to the innermost layer (at the interior wall surface).

Material
""""""""

Material composing the layer.

=============   ==================
**Required:**   Yes
**Type:**       Material reference
=============   ==================

Thickness
"""""""""

Thickness of the layer.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Interior Emissivity
^^^^^^^^^^^^^^^^^^^

Interior emissivity of the wall used for interior long-wave radiation calculations.

=============   =============
**Required:**   No
**Type:**       Numeric
**Units:**      dimensionless
**Default:**    0.8
=============   =============

Exterior Emissivity
^^^^^^^^^^^^^^^^^^^

Exterior emissivity of the wall used for exterior long-wave radiation calculations.

=============   =============
**Required:**   No
**Type:**       Numeric
**Units:**      dimensionless
**Default:**    0.8
=============   =============

Exterior Absorptivity
^^^^^^^^^^^^^^^^^^^^^

Exterior absorptivity of the wall used for calculating absorbed solar radiation.

=============   =============
**Required:**   No
**Type:**       Numeric
**Units:**      dimensionless
**Default:**    0.8
=============   =============

Interior Horizontal Insulation
------------------------------

This defines the position, dimensions, and material of interior horizontal insulation. Interior horizontal insulation begins at the wall’s interior surface and extends inward and downward to a user-specified width and thickness at a user-specified depth at or below the `Foundation Depth`_.

**Example:**

.. code-block:: yaml

  Interior Horizontal Insulation:
    Material: XPS # Material reference
    Thickness: 0.0508 # [m]
    Depth: 0.2032  # [m]
    Width: 0.4064 # [m]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============


Material
^^^^^^^^

Insulation material reference.

=============   ==================
**Required:**   Yes
**Type:**       Material reference
=============   ==================

Thickness
^^^^^^^^^

Thickness of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Depth
^^^^^

Depth of the insulation measured from the wall top to the top of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Width
^^^^^

Width of the insulation extending from the interior wall surface.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Interior Vertical Insulation
----------------------------

This defines the position, dimensions, and material of interior vertical insulation. Interior vertical insulation begins at the wall top and extends downward and inward to a user-specified depth and thickness. The depth can be specified to model partial interior wall insulation.

**Example:**

.. code-block:: yaml

  Interior Vertical Insulation:
    Material: XPS # Material reference
    Thickness: 0.0508 # [m]
    Depth: 0.6096  # [m]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============

Material
^^^^^^^^

Insulation material reference.

=============   ==================
**Required:**   Yes
**Type:**       Material reference
=============   ==================

Thickness
^^^^^^^^^

Thickness of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Depth
^^^^^

Depth of the insulation measured from the wall top to the bottom of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Exterior Horizontal Insulation
------------------------------

This defines the position, dimensions, and material of exterior horizontal insulation. Exterior horizontal insulation begins at the wall’s exterior surface and extends outward and downward to a user-specified width and thickness at a user-specified depth at or below the grade level.

**Example:**

.. code-block:: yaml

  Exterior Horizontal Insulation:
    Material: XPS # Material reference
    Thickness: 0.0508 # [m]
    Depth: 0.3048  # [m]
    Width: 0.6096 # [m]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============


Material
^^^^^^^^

Insulation material reference.

=============   ==================
**Required:**   Yes
**Type:**       Material reference
=============   ==================

Thickness
^^^^^^^^^

Thickness of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Depth
^^^^^

Depth of the insulation measured from the wall top to the top of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Width
^^^^^

Width of the insulation extending from the interior wall surface.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======


Exterior Vertical Insulation
----------------------------

This defines the position, dimensions, and material of exterior vertical insulation. Exterior vertical insulation begins at the wall top and extends downward and outward to a user-specified depth and thickness.

**Example:**

.. code-block:: yaml

  Exterior Vertical Insulation:
    Material: XPS # Material reference
    Thickness: 0.0508 # [m]
    Depth: 2.0  # [m]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============

Material
^^^^^^^^

Insulation material reference.

=============   ==================
**Required:**   Yes
**Type:**       Material reference
=============   ==================

Thickness
^^^^^^^^^

Thickness of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

Depth
^^^^^

Depth of the insulation measured from the wall top to the bottom of the insulation.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      m
=============   =======

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

If `Indoor Air Temperature Method`_ is ``FILE`` the indoor dry-bulb temperature will be set using hourly values defined in a comma separted value (CSV) file.

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

=============   ====
**Required:**   Yes
**Type:**       Path
=============   ====

Index
^^^^^

A list of two values corresponding to the row and column where the hourly data begins in the file. A value of ``[0, 0]`` starts at the first row and first column. A value of ``[0,1]`` starts at the first row and second column.

=============   ====================
**Required:**   Yes
**Type:**       List [2] of integers
=============   ====================

Indoor Air Temperature
----------------------

If `Indoor Air Temperature Method`_ is ``CONSTANT`` the indoor dry-bulb temperature will be set using this value. If `Indoor Air Temperature Method`_ is ``FILE``, then this is not required.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

Outdoor Temperature Method
--------------------------

Allows the user to choose between having a constant outdoor temperature for the duration of the simulaiton or to reference temperatures from the weather file.

=============   ================================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``WEATHER-FILE`` or ``CONSTANT``
**Default:**    ``WEATHER-FILE``
=============   ================================

Outdoor Dry-Bulb Temperature
----------------------------

If `Outdoor Temperature Method`_ is ``CONSTANT`` the outdoor dry-bulb temperature will be set using this value. If `Outdoor Temperature Method`_ is ``WEATHER-FILE``, then this is not required.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

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

Soil Absorptivity
-----------------

Solar absorptivity of the soil or grade surface.

=============   =============
**Required:**   No
**Type:**       Numeric
**Units:**      dimensionless
**Default:**    0.8
=============   =============

Soil Emissivity
---------------

Long-wave emissivity of the soil or grade surface.

=============   =============
**Required:**   No
**Type:**       Numeric
**Units:**      dimensionless
**Default:**    0.8
=============   =============

Surface Roughness
-----------------

Represents the relief of the surface. This value is used to calculate forced convection and the wind speed near the grade surface. Roughness values in the table below are converted from the more qualitative rougness values used in DOE-2 and EnergyPlus. Estimates for soil, gravel, and grass are also shown.

===============  =============
Example Surface  Roughness [m]
===============  =============
Glass                0.0000
Smooth Plaster       0.0044
Clear Pine           0.0052
Concrete             0.0208
Brick                0.0268
Stucco               0.0468
Soil                 0.0500
Gravel               0.1200
Grass                0.3000
===============  =============

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      m
**Default:**    0.3
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

Wall Top Temperature Difference
-------------------------------

If `Wall Top Boundary Condition`_ is ``LINEAR-DT``, then this value specifies the change in temperature across the wall. This is used only to represent the constraints of the IEA BESTEST analytical solution in case GC10a. The actual temperatures are determined based on the values of  `Indoor Air Temperature`_ and `Outdoor Dry-Bulb Temperature`_.

=============   =======
**Required:**   Depends
**Type:**       Numeric
**Units:**      K
=============   =======

Orientation
-----------

Defines the orientation of the building clockwise relative to North (East = :math:`\pi/2`, South = :math:`\pi`, West = :math:`3\pi/2`). This is used to calculate the solar incidence and wind direction relative to exterior vertical foundation surfaces.

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      radians
**Default:**    0.0
=============   =======


Number of Dimensions
--------------------

Coordinate System
-----------------

Two-Dimensional Approximation
-----------------------------

Length 1
--------

Length 2
--------

Use Symmetry
------------

Perimeter Surface Width
-----------------------

Mesh
----

Minimum Cell Dimension
^^^^^^^^^^^^^^^^^^^^^^

Maximum Near-Field Growth Coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Maximum Deep-Field Growth Coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Maximum Interior-Field Growth Coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Maximum Far-Field Growth Coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Numerical Scheme
----------------

f-ADI
-----

Solver
------

Preconditioner
--------------

Maximum Iterations
------------------

Tolerance
---------

Initialization Method
---------------------

Initial Temperature
-------------------

Accelerated Initialization Timestep
-----------------------------------

Days

Number of Accelerated Initialization Timesteps
----------------------------------------------

Number of Warmup Days in Initialization
---------------------------------------

Output Report
-------------

Reports
^^^^^^^

Minimum Reporting Frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Output Snapshots
----------------

Name
^^^^

Size
^^^^

Frequency
^^^^^^^^^

Start Date
^^^^^^^^^^

End Date
^^^^^^^^

X Range
^^^^^^^

Y Range
^^^^^^^

Z Range
^^^^^^^

File Format
^^^^^^^^^^^

Unit System
^^^^^^^^^^^

Plot Type
^^^^^^^^^

Flux Direction
^^^^^^^^^^^^^^

Output Range
^^^^^^^^^^^^



Color Scheme
^^^^^^^^^^^^

Grid
^^^^

Axes
^^^^

Timestamp
^^^^^^^^^

Gradients
^^^^^^^^^

Contours
^^^^^^^^

Contour Labels
^^^^^^^^^^^^^^

Number of Contours
^^^^^^^^^^^^^^^^^^
