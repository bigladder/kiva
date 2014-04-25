#!/bin/bash
params compose ../../../templates/template.params -f "../../test.txt;../case.txt;soln.txt" -o instance.yaml
../../../build/kiva instance.yaml ../*.epw Timeseries.csv
