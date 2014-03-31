#!/bin/bash
params compose ../../../../templates/template.params -f params.txt -o instance.yaml
../../../../build/kiva instance.yaml ../*.epw Timeseries.csv
