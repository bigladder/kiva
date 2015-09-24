Materials
=========

Materials define the intrensic thermophysical properties of the solid components of the computational domain (eg: surrounding soil, structural concrete, insulation). All materials must fall within a single `Materials:` group.

**Example:**

.. code-block:: yaml

    Materials:
      Soil:
        Conductivity: 0.864 # [W/m-K]
        Density: 1510.0  # [kg/m3]
        Specific Heat: 1260.0  # [J/kg-K]
      Concrete:
        Conductivity: 1.98  # [W/m-K]
        Density: 1900.0  # [kg/m3]
        Specific Heat: 665.0  # [J/kg-K]
      XPS:
        Conductivity: 0.029  # [W/m-K]
        Density: 28.0  # [kg/m3]
        Specific Heat: 1450.0  # [J/kg-K]

Each instance of a material begins with a descriptive name (eg: Soil, Concrete or XPS) and contains the following attributes:

Conductivity
------------

Thermal conductivity of the material.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      W/m-K
=============   =======

Density
-------

Density of the material.

=============   ==============
**Required:**   Yes
**Type:**       Numeric
**Units:**      kg/m\ :sup:`3`
=============   ==============

Specific Heat
-------------

Specific heat of the material.

=============   =======
**Required:**   Yes
**Type:**       Numeric
**Units:**      J/kg-K
=============   =======
