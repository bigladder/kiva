import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import xlrd
import datetime as dt
import csv

from dateutil.parser import parse

def getLastValue(case, solution):
    with open('../'+case+'/'+solution+'/Timeseries.csv') as f:
        lines = f.readlines()
        reader = csv.reader([lines[-1]])
        for row in reader:
            val = row[1]

    return float(val)

print "Begin generating plots..."

sns.set_style("white")

# GC10a

# Absolute values

fig, ax = plt.subplots()

refs = ['Analytical','TRNSYS', 'FLUENT', 'MATLAB']

blue_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[0]],5)
green_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[1]],5)
red_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[2]],5)
purple_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[3]],5)
yellow_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[4]],5)

colors = [blue_palette[4],
          green_palette[4],
          green_palette[3],
          green_palette[2]]

gc10a_values = []

print "Read reference results..."
workbook = xlrd.open_workbook('../doc/GC-InDepth-Results.XLS')

for name in refs:
    worksheet = workbook.sheet_by_name(name)
    gc10a_values.append(worksheet.cell_value(57,4))

names = refs

gc10a_results = {}


# Analytical
'''
val = gc10a_values[0]

dates = [parse('2096-Dec-31 00:00:00'), parse('2097-Dec-31 23:00:00')]
ts = pd.Series([val,val], index=dates)
gc10a_results['Analytical'] = pd.DataFrame({'W': ts})
'''
# SS
print "Read SS results..."

df = pd.read_csv('../GC10a/SteadyState/Timeseries.csv',
                                   header=0,
                                   names=['time','W'],
                                   parse_dates=True,
                                   index_col=0)
names.append('Kiva: SS')

val = df['W'][len(df.index)-1]
gc10a_values.append(val)

'''
dates = [parse('2096-Dec-31 00:00:00'), parse('2097-Dec-31 23:00:00')]
ts = pd.Series([val,val], index=dates)
gc10a_results['SS'] = pd.DataFrame({'W': ts})
'''

colors.append(purple_palette[4])

# ADE
print "Read ADE results..."

gc10a_results['ADE'] = pd.read_csv('../GC10a/ADE/Timeseries.csv',
                                   header=0,
                                   names=['time','W'],
                                   parse_dates=True,
                                   index_col=0)
df = gc10a_results['ADE']
names.append('Kiva: ADE')
gc10a_values.append(df['W'][len(df.index)-1])
colors.append(red_palette[4])

# ADI
print "Read ADI results..."

gc10a_results['ADI'] = pd.read_csv('../GC10a/ADI/Timeseries.csv',
                                   header=0,
                                   names=['time','W'],
                                   parse_dates=True,
                                   index_col=0)
df = gc10a_results['ADI']
names.append('Kiva: ADI')
gc10a_values.append(df['W'][len(df.index)-1])
colors.append(red_palette[3])

# Implicit
print "Read Implicit results..."

'''
gc10a_results['Implicit'] = pd.read_csv('../GC10a/Implicit/Timeseries.csv',
                                   header=0,
                                   names=['time','W'],
                                   parse_dates=True,
                                   index_col=0)
df = gc10a_results['Implicit']
'''

names.append('Kiva: Implicit')
gc10a_values.append(df['W'][len(df.index)-1])
colors.append(yellow_palette[4])

# Explicit
print "Read Explicit results..."

'''
gc10a_results['Explicit'] = pd.read_csv('../GC10a/Explicit/Timeseries.csv',
                                   header=0,
                                   names=['time','W'],
                                   parse_dates=True,
                                   index_col=0)
df = gc10a_results['Explicit']
'''

names.append('Kiva: Explicit')
gc10a_values.append(getLastValue('GC10a','Explicit'))
colors.append(yellow_palette[3])

#print gc10a_values

# Bar Chart
print "Create figure..."

ind = range(1,len(names)+1)
width = 1

data = ax.bar(ind, gc10a_values, width, color=colors)

'''
# Steady-State time series chart

df = pd.concat([gc10a_results['Analytical'].resample('H', fill_method='pad'),
                gc10a_results['SS'].resample('H', fill_method='pad'),
                gc10a_results['ADE'],
                gc10a_results['ADI'],
                gc10a_results['Implicit'],
                gc10a_results['Explicit'].resample('H', how='mean')], 
               axis=1, keys=['Analytical','Steady State','ADE','ADI','Implicit','Explicit'])

df.plot()
'''

plt.show()

print "Done."



if __name__ == '__main__':
    pass