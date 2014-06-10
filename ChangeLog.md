0.2.0 Released ? ??? 2014
-------------------------

### Features:
* Change to LIS (Library of Iterative Solvers) for linear algebra solutions. User can
  select solver options.
* Add ability to model partial wall insulation for basements.
* Add ability to simulate shading using graphics hardware acceleration.
* New initialization scheme separates initialization into initial conditions, implicit
  acceleration and warm-up periods. All three may be used for a single initialization.
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
* Fix interpretation of weather file data.
* Implicit acceleration is now less memory-intensive.
* Use Tridiagonal Matrix Algorithm for ADI calculations (dramatically improves runtime).
* Fix windward convection algorithm typo.
* Fix plot update for plots starting mid-year.

0.1.1 Released 9 December 2013
------------------------------
* Make output file a command line argument.

0.1.0 Released 6 December 2013
------------------------------
* Initial release.