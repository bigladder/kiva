#!/bin/bash

echo "Kiva(TM) is a command line program. To run it, you must type on the command prompt."
echo ""
echo "You can see how it's used by typing:"
echo ""
echo "./kiva --help"
echo ""
read -n1 -r -p "Press any key to see the output . . ."
echo ""

./kiva --help

read -n1 -r -p "Press any key to continue . . ."
echo ""
echo "As you can see, Kiva takes three inputs:"
echo "  1. A YAML input file describing the foundation and simulation settings"
echo "  2. An EnergyPlus weather file (EPW) representing your location"
echo "  3. The name of your output file (CSV)"
echo ""
echo "For example, to simulate the slab example in Golden, Colorado, type:"
echo ""
echo "./kiva examples/slab.yaml weather/USA_CO_Golden.epw kiva-output.csv"
echo ""
read -n1 -r -p "Press any key to run this example in Kiva . . ."
echo ""

./kiva examples/slab.yaml weather/USA_CO_Golden.epw kiva-output.csv

echo ""
echo "You now have a CSV output file called 'kiva-output.csv' that you can open"
echo "in a spreadsheet!"
echo ""
echo "For more details see the user's guide: http://kiva.readthedocs.io/"
echo ""
read -n1 -r -p "Press any key to leave the demo . . ."
