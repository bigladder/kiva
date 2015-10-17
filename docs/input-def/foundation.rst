.. _foundation:

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

.. _polygon:

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

Orientation
-----------

Defines the orientation of the building clockwise relative to North (East = :math:`\pi/2`, South = :math:`\pi`, West = :math:`3\pi/2`). This is used to calculate the solar incidence and wind direction relative to exterior vertical foundation surfaces.

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      radians
**Default:**    0.0
=============   =======

Perimeter Surface Width
-----------------------

This value is used to define a portion of the slab's perimeter separately from the slab core. This will affect the meshing of the slab, but is intended primarily for separate output reporting for each region.

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      m
**Default:**    0.0
=============   =======

Output Report
-------------

The output report defines what variables are written to the CSV output file and how often they are written.

**Example:**

.. code-block:: yaml

  Output Report:
    Minimum Reporting Frequency: 60 # [min]
    Reports:
      - 0 # Slab Core Average Heat Flux [W/m2]
      - 1 # Slab Core Average Temperature [K]
      - 2 # Slab Core Average Effective Temperature [C]
      - 3 # Slab Core Total Heat Transfer Rate [W]
      - 4 # Slab Perimeter Average Heat Flux [W/m2]
      - 5 # Slab Perimeter Average Temperature [K]
      - 6 # Slab Perimeter Average Effective Temperature [C]
      - 7 # Slab Perimeter Total Heat Transfer Rate [W]
      - 8 # Slab Average Heat Flux [W/m2]
      - 9 # Slab Average Temperature [K]
      - 10 # Slab Total Heat Transfer Rate [W]
      - 11 # Wall Average Heat Flux [W/m2]
      - 12 # Wall Average Temperature [K]
      - 13 # Wall Average Effective Temperature [C]
      - 14 # Wall Total Heat Transfer Rate [W]
      - 15 # Foundation Average Heat Flux [W/m2]
      - 16 # Foundation Average Temperature [K]
      - 17 # Foundation Total Heat Transfer Rate [W]

=============   ===============
**Required:**   No
**Type:**       Compound object
=============   ===============

Minimum Reporting Frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Kiva can provide timeseries output at the same interval as the timestep. This input will override to delay output and write it a a lower frequency. This helps to reduce the output size when running at very small timesteps.

=============   =======
**Required:**   No
**Type:**       Integer
**Units:**      min
**Default:**    60
=============   =======

Reports
^^^^^^^

This is a list of report ID numbers that Kiva will write to the CSV output file. The IDs and there corresponding output are listed in the table below:

====  ===========================================   ===============
ID    Output Variable                               Units
====  ===========================================   ===============
0     Slab Core Average Heat Flux                   W/m\ :sup:`2`
1     Slab Core Average Temperature                 K
2     Slab Core Average Effective Temperature       :math:`^\circ`C
3     Slab Core Total Heat Transfer Rate            W
4     Slab Perimeter Average Heat Flux              W/m\ :sup:`2`
5     Slab Perimeter Average Temperature            K
6     Slab Perimeter Average Effective Temperature  :math:`^\circ`C
7     Slab Perimeter Total Heat Transfer Rate       W
8     Slab Average Heat Flux                        W/m\ :sup:`2`
9     Slab Average Temperature                      K
10    Slab Total Heat Transfer Rate                 W
11    Wall Average Heat Flux                        W/m\ :sup:`2`
12    Wall Average Temperature                      K
13    Wall Average Effective Temperature            :math:`^\circ`C
14    Wall Total Heat Transfer Rate                 W
15    Foundation Average Heat Flux                  W/m\ :sup:`2`
16    Foundation Average Temperature                K
17    Foundation Total Heat Transfer Rate           W
====  ============================================  ===============

When `Perimeter Surface Width`_ is not specified, the entire slab is considered to be "Core".

"Effective Temperature" is used for preprocessed ground temperatures in whole-building simulation engines. These values represent the effective temperature on the ground's side of the slab core, slab perimeter, or wall layers. When used in a whole-building simulation, the construction in the whole-building model should be the same as the layers defined for the respective surface in Kiva (ignoring any insulation objects).

=============   ====================
**Required:**   No
**Type:**       List [N] of integers
**Default:**    No reports
=============   ====================

Output Snapshots
----------------

Output snapshots are used to graphically visualize domain temperatures and/or heat fluxes. Each series of snapshots is part of a list within the `Output Snapshots`_ object. A series consists of potentially many snapshots taken of a slice of the domain at a user-specified frequency between a start and end date.

.. figure:: ../images/snapshot-profile.png

  Example profile snapshot

.. figure:: ../images/snapshot-plan.png

  Example plan snapshot

**Example:**

.. code-block:: yaml

  Output Snapshots:
    -
     Directory: Output/Profile
     Size: 800
     Frequency: 1
     Start Date: 2015-Dec-21
     End Date: 2015-Dec-21
     X Range: [0, 30]
     Z Range: [-30, 0.3048]



=============   ============================
**Required:**   No
**Type:**       List [N] of compound objects
=============   ============================

Directory
^^^^^^^^^

Directory where snapshots are created. An ordered file name, ``XXXX.png``, identifies each snapshot within a series. For example, the 134th snapshot in a series with a directory name of ``Profile`` will be created as ``Profile/0134.png``.

=============   ==============
**Required:**   Yes
**Type:**       Directory Path
=============   ==============

Size
^^^^

The size in pixels of each snapshot file. Outputs are all generated as square images.

=============   =======
**Required:**   No
**Type:**       Integer
**Units:**      pixels
**Default:**    800
=============   =======

Frequency
^^^^^^^^^

The frequency, in hours, at which new snapshots are taken. The default is 36 hours so that the snapshots capture both nighttime and daytime output.

=============   =======
**Required:**   No
**Type:**       Integer
**Units:**      hours
**Default:**    36
=============   =======

Start Date
^^^^^^^^^^

Specifies the start date of the snapshots. Snapshots begin at 12:00am of this day. This is specified as a date string (e.g., YYYY-Mon-DD, YYYY/MM/DD).

=============   =====================
**Required:**   No
**Type:**       Date string
**Default:**    Simulation start date
=============   =====================

End Date
^^^^^^^^

Specifies the end date of the snapshots. Snapshots end before 12:00am of the following day. This is specified as a date string (e.g., YYYY-Mon-DD, YYYY/MM/DD).

=============   ===================
**Required:**   No
**Type:**       Date string
**Default:**    Simulation end date
=============   ===================

X Range
^^^^^^^

Defines the range the domain captured in the snapshot in the "X"-direction (``[Xmin, Xmax]``). By default the `X Range`_ will show the entire extents of the "X" direction, and may not show the detail where heat is flowing near the foundaiton. For three-dimensional solutions, a slice along a plane in the "X"-direction can be specified by giving both ``Xmin`` and ``Xmax`` the same value.

The snapshot will round the range to the next cell division.

=============   =========================
**Required:**   No
**Type:**       List [2] of numerics
**Units:**      m
**Default:**    "X" extents of the domain
=============   =========================

Y Range
^^^^^^^

Defines the range the domain captured in the snapshot in the "Y"-direction (``[Ymin, Ymax]``). By default the `Y Range`_ will show the entire extents of the "Y" direction, and may not show the detail where heat is flowing near the foundaiton. For three-dimensional solutions, a slice along a plane in the "Y"-direction can be specified by giving both ``Ymin`` and ``Ymax`` the same value. For two-dimensional simulations this should not be included.

The snapshot will round the range to the next cell division.

=============   =========================
**Required:**   No
**Type:**       List [2] of numerics
**Units:**      m
**Default:**    "Y" extents of the domain
=============   =========================

Z Range
^^^^^^^

Defines the range the domain captured in the snapshot in the "Z"-direction (``[Zmin, Zmax]``). By default the `Z Range`_ will show the entire extents of the "Z" direction, and may not show the detail where heat is flowing near the foundaiton. For three-dimensional solutions, a slice along a plane in the "Z"-direction can be specified by giving both ``Zmin`` and ``Zmax`` the same value.

The snapshot will round the range to the next cell division.

=============   =========================
**Required:**   No
**Type:**       List [2] of numerics
**Units:**      m
**Default:**    "Z" extents of the domain
=============   =========================

Plot Type
^^^^^^^^^

Defines the type of output plotted. Options are ``TEMPERATURE`` and ``HEAT-FLUX``. For ``HEAT-FLUX``, the user may also specify a `Flux Direction`_ for output.

=============   ================================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``TEMPERATURE`` or ``HEAT-FLUX``
**Default:**    ``TEMPERATURE``
=============   ================================

Flux Direction
^^^^^^^^^^^^^^

When `Plot Type`_ is ``HEAT-FLUX``, the snapshots show the magnitude of heat flux throughout the domain. This input allows the user to specify whether they want to display the overall magnitude, ``MAG``, or the magnitude in a given direciton, ``X``, ``Y``, or ``Z``.

=============   =======================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``MAG``, ``X``, ``Y``, or ``Z``
**Default:**    ``MAG``
=============   =======================

Unit System
^^^^^^^^^^^

Defines the units used in the output snapshots. Options are ``IP`` (Inch-Pound), and ``SI`` (International System). Keep in mind that regardless of this value, all other inputs are still defined in the SI unit system.

=============   ================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``IP`` or ``SI``
**Default:**    ``SI``
=============   ================

Output Range
^^^^^^^^^^^^

Specifies the range of output shown in the snapshots. The units of the range depend on the value of `Plot Type`_ and `Unit System`_.

=============   ==================================
**Required:**   No
**Type:**       List [2] of numerics
**Units:**      Depends
**Default:**    "Z" extents of the domain
=============   ==================================

Color Scheme
^^^^^^^^^^^^

Specifies the color scheme used within the `Output Range`_. Options are:

- ``CMR``, best color scheme where colors progress in brightness with magnitude (prints in black-and-white),
- ``JET``, like a rainbow(!), but doesn't print well,
- ``NONE``, do not show any output. This can be used to illustrate meshing independent of results.

=============   =============================
**Required:**   No
**Type:**       Enumeration
**Values:**     ``CMR``, ``JET``, or ``NONE``
**Default:**    ``CMR``
=============   =============================

Mesh
^^^^

Enables the display of the mesh (discretized cells).

=============   =======
**Required:**   No
**Type:**       Boolean
**Default:**    False
=============   =======

Axes
^^^^

Enables the display of the spatial axes, and the colorbar.

=============   =======
**Required:**   No
**Type:**       Boolean
**Default:**    True
=============   =======

Timestamp
^^^^^^^^^

Enables the display of the timestamp.

=============   =======
**Required:**   No
**Type:**       Boolean
**Default:**    True
=============   =======

Gradients
^^^^^^^^^

Enables the display of gradients.

=============   =======
**Required:**   No
**Type:**       Boolean
**Default:**    False
=============   =======

Contours
^^^^^^^^

Enables the display of contours.

=============   =======
**Required:**   No
**Type:**       Boolean
**Default:**    True
=============   =======

Contour Labels
^^^^^^^^^^^^^^

Enables the display of contour labels.

=============   =======
**Required:**   No
**Type:**       Boolean
**Default:**    False
=============   =======

Number of Contours
^^^^^^^^^^^^^^^^^^

Specifies the number of countours to generate between the values specified in `Output Range`_.

=============   =======
**Required:**   No
**Type:**       Integer
**Default:**    13
=============   =======
