@echo off
echo Kiva(TM) is a command line program. To run it, you must type on the command
echo prompt.
echo.
echo You can see how it's used by typing:
echo.
echo kiva.exe --help
echo.
echo Press any key to see the output . . .
pause >nul
echo.

kiva.exe --help

pause
echo.
echo As you can see, Kiva takes three inputs:
echo   1. A YAML input file describing the foundation and simulation settings
echo   2. An EnergyPlus weather file (EPW) representing your location
echo   3. The name of your output file (CSV)
echo.
echo For example, to simulate the slab example in Golden, Colorado, type:
echo.
echo kiva.exe examples/slab.yaml weather/USA_CO_Golden.epw kiva-output.csv
echo.
echo Press any key to run this example in Kiva . . .
pause >nul
echo.

kiva.exe examples/slab.yaml weather/USA_CO_Golden.epw kiva-output.csv

echo.
echo You now have a CSV output file called 'kiva-output.csv' that you can open
echo in a spreadsheet!
echo.
echo For more details see the user's guide: http://kiva.readthedocs.io/
echo.
echo Press any key to leave the demo . . .
pause >nul
