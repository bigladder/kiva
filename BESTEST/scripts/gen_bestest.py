print "Import libraries..."

print "Import matplotlib..."

import matplotlib.pyplot as plt

print "Import seaborn..."

import seaborn as sns

print "Import xlrd..."

import xlrd

print "Import csv..."

import csv

print "Read Excel results..."


workbook = xlrd.open_workbook('../doc/GC-InDepth-Results.XLS')

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

class Solution:
    def __init__(self):
        self.name = ""
        self.color = []
        self.hatch = ""
        self.ID = ""
        self.ref = False
        
print "Read results..."    

sns.set_style("whitegrid", {'axes.grid': False})
blue_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[0]],5)
green_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[1]],5)
red_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[2]],5)
purple_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[3]],5)
yellow_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[4]],5)

'''
colors = [blue_palette[4],
          green_palette[4],
          green_palette[3],
          green_palette[2],
          purple_palette[4],
          red_palette[4],
          red_palette[3],
          yellow_palette[4],
          yellow_palette[3],
          yellow_palette[2]]
'''

colors = ['m',
          '#1BC3F9',
          '#1BC3F9',
          '#1BC3F9',
          'w',
          'w',
          'w',
          'w',
          'w',
          'w']

hatches = ['',
           'ooo',
           'xxxx',
           '.....',
           '////',
           '\\\\\\\\',
           'xxxx',
           '---',
           '------',
           '+++']

solutions = []

fig = plt.figure()
ax = plt.subplot(111)

j = 0

analytical_solution = Solution()
analytical_solution.color = colors[j]
analytical_solution.hatch = hatches[j]
analytical_solution.ref = True
analytical_solution.name = "Analytical"
analytical_solution.ID = "Analytical"
solutions.append(analytical_solution)
j += 1

trnsys_solution = Solution()
trnsys_solution.color = colors[j]
trnsys_solution.hatch = hatches[j]
trnsys_solution.ref = True
trnsys_solution.name = "TRNSYS"
trnsys_solution.ID = "TRNSYS"
solutions.append(trnsys_solution)
j += 1

fluent_solution = Solution()
fluent_solution.color = colors[j]
fluent_solution.hatch = hatches[j]
fluent_solution.ref = True
fluent_solution.name = "FLUENT"
fluent_solution.ID = "FLUENT"
solutions.append(fluent_solution)
j += 1

matlab_solution = Solution()
matlab_solution.color = colors[j]
matlab_solution.hatch = hatches[j]
matlab_solution.ref = True
matlab_solution.name = "MATLAB"
matlab_solution.ID = "MATLAB"
solutions.append(matlab_solution)
j += 1

ss_solution = Solution()
ss_solution.color = colors[j]
ss_solution.hatch = hatches[j]
ss_solution.ref = False
ss_solution.name = "Kiva: Steady-State"
ss_solution.ID = "SteadyState"
solutions.append(ss_solution)
j += 1

ade_solution = Solution()
ade_solution.color = colors[j]
ade_solution.hatch = hatches[j]
ade_solution.ref = False
ade_solution.name = "Kiva: ADE"
ade_solution.ID = "ADE"
solutions.append(ade_solution)
j += 1

adi_solution = Solution()
adi_solution.color = colors[j]
adi_solution.hatch = hatches[j]
adi_solution.ref = False
adi_solution.name = "Kiva: ADI"
adi_solution.ID = "ADI"
solutions.append(adi_solution)
j += 1

implicit_solution = Solution()
implicit_solution.color = colors[j]
implicit_solution.hatch = hatches[j]
implicit_solution.ref = False
implicit_solution.name = "Kiva: Implicit"
implicit_solution.ID = "Implicit"
solutions.append(implicit_solution)
j += 1

explicit_solution = Solution()
explicit_solution.color = colors[j]
explicit_solution.hatch = hatches[j]
explicit_solution.ref = False
explicit_solution.name = "Kiva: Explicit"
explicit_solution.ID = "Explicit"
solutions.append(explicit_solution)
j += 1

cn_solution = Solution()
cn_solution.color = colors[j]
cn_solution.hatch = hatches[j]
cn_solution.ref = False
cn_solution.name = "Kiva: Crank-Nicolson"
cn_solution.ID = "CrankNicolson"
solutions.append(cn_solution)
j += 1

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

for soln in solutions:
    names.append(soln.name)
    colors.append(soln.color)
    hatches.append(soln.hatch)

ind = []

data = []

width = 1

ticks = []
i = 1

# Bar Chart
print "Create figure..."

for case in cases:
    ticks.append(i+5)
    for soln in solutions:
        values.append(getValue(case,soln))
        data.append(ax.bar(i, getValue(case,soln), width, color=soln.color, hatch=soln.hatch))
        ind.append(i)
        i+=1
    i+=1
    


box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
ax.set_autoscalex_on(False)
ax.set_xlim([0,i])
#ax.set_title('IEA BESTEST Ground Coupling: In-Depth Floor Slab\nSteady-State Floor Conduction')
ax.set_xticks(ticks)
ax.set_xticklabels(cases)
ax.yaxis.grid()

legend = ax.legend(data[:10], names[:10], loc='upper center', ncol=5, bbox_to_anchor=(0.5,-0.05),
                   fancybox=True)


#plt.show()
plt.savefig('../figures/BESTEST_SS.pdf')

print "Done."



if __name__ == '__main__':
    pass