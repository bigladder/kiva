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

Layers are specified as a list of material references, and thickness values from the outtermost layer to the innermost layer (at the floor surface).

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

=============   =======
**Required:**   No
**Type:**       Numeric
**Units:**      m
**Default:**    0.8
=============   =======

Wall
----

This defines the costruction of the foundation wall. This is not required. If there is no wall defined for a given foundation, then the wall is exposed soil.

**Example:**

.. code-block:: yaml

  Wall:
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

defines the height of the wall top relative to the grade (z = 0).

Layers
^^^^^^

Layers are specified as a list of material references, and thickness values from the outtermost layer to the innermost layer (at the interior wall surface).

Material
""""""""

Thickness
"""""""""

Interior Emissivity
^^^^^^^^^^^^^^^^^^^

Exterior Emissivity
^^^^^^^^^^^^^^^^^^^

Exterior Absorptivity
^^^^^^^^^^^^^^^^^^^^^

Interior Horizontal Insulation
------------------------------

Material
^^^^^^^^

Thickness
^^^^^^^^^

Depth
^^^^^

Width
^^^^^

Interior Vertical Insulation
----------------------------

Material
^^^^^^^^

Thickness
^^^^^^^^^

Depth
^^^^^

Exterior Horizontal Insulation
------------------------------

Material
^^^^^^^^

Thickness
^^^^^^^^^

Depth
^^^^^

Width
^^^^^

Exterior Vertical Insulation
----------------------------

Material
^^^^^^^^

Thickness
^^^^^^^^^

Depth
^^^^^

Indoor Air Temperature Method
-----------------------------

Indoor Air Temperature File
---------------------------

Name
^^^^

Index
^^^^^

List of two values, row/column.

Indoor Air Temperature
----------------------

Outdoor Temperature Method
--------------------------

Outdoor Dry-Bulb Temperature
----------------------------


Far-Field Width
---------------

Deep-Ground Depth
-----------------

Deep-Ground Boundary Condition
------------------------------

Convection Calculation Method
-----------------------------

Interior Convection Coefficient
-------------------------------

Exterior Convection Coefficient
-------------------------------

Soil Absorptivity
-----------------

Soil Emissivity
---------------

Surface Roughness
-----------------

Vegetation Height
-----------------

Delta Local
-----------

Alpha Local
-----------

Wall Top Boundary Condition
---------------------------

Wall Top Temperature Difference
-------------------------------

Orientation
-----------


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
