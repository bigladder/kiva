0.4.0 Released ????
-------------------
### Fixes:
* Improve output reporting. 2D simulations now report heat transfer rates
  representative of the 3D inputs.
* Allow interior temperature file to be relative to input file (or in working
  directory).
* Create output snapshots directory in same directory as the simulation output
  CSV file.
* Small fixes to solar and convective boundary conditions.
* Better error handling.

### Features:
* Ability to define any block of material within the domain.
* Add inputs to define the foundation footing.
* Add ability to define exposed foundation perimeters (including unexposed/core
  foundations).
* Make horizontal insulation depths relative to the slab and grade surfaces.
* Replace foundation wall height ("Height") with depth relative to bottom
  of the slab ("Depth Below Slab").
* Move boundary, initialization, numerical settings, and output related inputs
  to a higher level.
* Switch to simpler linear solver ([Eigen](http://eigen.tuxfamily.org/)).
  Removes inputs for solver and preconditioner.
* Separate core functionality into a linked libraries.
* Use continuous integration for building and testing.

0.3.1 Released 16 October 2015
------------------------------
### Fixes:
* Make indoor/outdoor air temperature fields more consistent.
* Cross-platform CSV line ending support.

0.3.0 Released 28 September 2015
--------------------------------
### Features:
* Better naming conventions for input fields.
* More informed default values.
* User documentation.
* Example files.

0.2.1 Released 16 March 2015
----------------------------
### Features:
* Introduces new two-dimensional approximation methods.
* Added ability to plot heat flux in addition to temperature.
* Added the ability to read indoor air temperatures from a CSV file.
* Microsecond resolution for elapsed time listed in stdout.
* More output variables:
  * Slab Average Heat Flux [W/m2]
  * Slab Average Temperature [K]
  * Foundation Average Heat Flux [W/m2]
  * Foundation Average Temperature [K]

### Fixes:
* Fixes for symmetric simulations.

0.2.0 Released 11 June 2014
---------------------------
### Features:
* Change to LIS (Library of Iterative Solvers) for linear algebra solutions.
  User can select solver options.
* Add ability to model partial wall insulation for basements.
* Add ability to simulate shading using graphics hardware acceleration.
* New initialization scheme separates initialization into initial conditions,
  implicit acceleration and warm-up periods. All three may be used for a single
  initialization.
* Allow user to specify the minimum output frequency (e.g., output every hour
  even if timestep is smaller).
* Generalize hourly output format.
* Allow user to specify exterior grade roughness (e.g., for grass or asphalt).

### Cosmetic Changes:
* Make CMake builds platform independent.
* Executable is now lower case: "kiva".
* Move all threading to OpenMP.
* Provide better status updates to user during initialization.
* Change plot outputs to IP units and add temperature units to color bar.
* Change wall dimension variable to "height above grade" (was "wall depth").
* Change "excavation depth" variable to "foundation depth".

### Fixes:
* Calculate incident solar on vertical surfaces (previously used horizontal
  value).
* Fix interpretation of weather file data.
* Implicit acceleration is now less memory-intensive.
* Use Tridiagonal Matrix Algorithm for ADI calculations (dramatically improves
  runtime).
* Fix problem where interior cells were not growing towards the center in 3D
  simulations.
* Fix windward convection algorithm typo.
* Fix block outlines in 3D slice plots.
* Fix plot update for plots starting mid-year.

0.1.1 Released 9 December 2013
------------------------------
* Make output file a command line argument.

0.1.0 Released 6 December 2013
------------------------------
* Initial release.
