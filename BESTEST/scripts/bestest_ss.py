print "Import libraries..."

from gen_bestest import *

def getValue(case, solution):
    if (solution.ref):
        return readXLS(case,solution.ID)
    else:
        return getLastValue(case, solution.ID)


def readXLS(case, ID):
    worksheet = workbook.sheet_by_name(ID)
    if (case == 'GC10a'):
        row = 57
    elif (case == 'GC30a'):
        row = 58
    elif (case == 'GC30b'):
        row = 59
    elif (case == 'GC30c'):
        row = 60
    elif (case == 'GC60b'):
        row = 61
    elif (case == 'GC65b'):
        row = 62
    
    if (ID == 'Analytical' and case != 'GC10a'):
        return 0.0
    else:    
        return worksheet.cell_value(row,4)
        
def getLastValue(case, ID):
    with open('../'+case+'/'+ID+'/Timeseries.csv') as f:
        lines = f.readlines()
        reader = csv.reader([lines[-1]])
        for row in reader:
            val = row[1]

    return float(val)

def fmt(x):
    return '{:,.0f}'.format(x)

print "Read results..."    

cases = ['GC10a',
         'GC30a',
         'GC30b',
         'GC30c',
         'GC60b',
         'GC65b']

values = []
names = []
colors = []
hatches = []
times = []

for soln in solutions:
    names.append(soln.name)
    colors.append(soln.color)
    hatches.append(soln.hatch)

data = []

time_data = []

width = 1

ticks = []
ticks2 = []
i = 1
i2 = 1

print "Create figures..."
file_name = 'bestest_ss'

# Figure style
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_context("paper", {'axes.labelsize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

fig = plt.figure()
ax = fig.add_subplot(111)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

for case in cases:
    ticks.append(i+5)
    ticks2.append(i2+3)
    for soln in solutions:
        value = getValue(case,soln)
        values.append(value)
        if value > 0:
            soln.values.append(value)
        else:
            soln.values.append('--')
        data.append(ax.bar(i, value, width, color=soln.color, hatch=soln.hatch))
        time = getTime(case,soln)
        soln.times.append(time)
        times.append(time)
        if (not soln.ref):
            time_data.append(ax2.bar(i2, time, width, color=soln.color, hatch=soln.hatch, log=True))
            i2 += 1
        i+=1
    i+=1
    i2+=1

# Bar Chart

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
ax.set_autoscalex_on(False)
ax.set_xlim([0,i])
#ax.set_title('IEA BESTEST Ground Coupling: In-Depth Floor Slab\nSteady-State Floor Conduction')
ax.set_xticks(ticks)
ax.set_xticklabels(cases)
ax.yaxis.grid()
ax.set_ylabel('Floor Heat Flow [W]')

legend = ax.legend(data[:10], names[:10], loc='upper center', ncol=5, bbox_to_anchor=(0.5,-0.05),
                   fancybox=True)



fig.savefig('../figures/' + file_name + '.pdf')
fig.savefig(output_dir + 'images/' + file_name + '.pdf')



# Time compare
box = ax2.get_position()
ax2.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

ax2.set_xlim([0,i2])
ax2.set_ylim([1,100000])
ax2.set_xticks(ticks2)
ax2.set_xticklabels(cases)
ax2.yaxis.grid()
ax2.set_ylabel('Simulation Time [s]')

legend2 = ax2.legend(time_data[:6], names[4:], loc='upper center', ncol=3, bbox_to_anchor=(0.5,-0.05),
                   fancybox=True)


fig2.savefig(output_dir + 'images/' + file_name + '_times.pdf')


# Create Table
print "Create table..."

d = {}
for soln in solutions:
    d[soln.name] = soln.values
    
df = pd.DataFrame(d, index=cases)
df = df[names].T

df.to_latex(output_dir + 'tables/' + file_name + '.tex',float_format=fmt)

d = {}
for soln in solutions:
    if (not soln.ref):
        d[soln.name] = soln.times
    
df = pd.DataFrame(d, index=cases)
df = df[names[4:]].T

df.to_latex(output_dir + 'tables/' + file_name + '_times.tex',float_format=fmt)




print "Done."



if __name__ == '__main__':
    pass