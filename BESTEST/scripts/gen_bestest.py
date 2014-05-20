import matplotlib.pyplot as plt
import seaborn as sns
import xlrd
import csv

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
        self.ID = ""
        self.ref = False
        
    

print "Begin generating plots..."

sns.set_style("white")
blue_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[0]],5)
green_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[1]],5)
red_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[2]],5)
purple_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[3]],5)
yellow_palette = sns.blend_palette(["ghostwhite",sns.color_palette("deep",6)[4]],5)


fig = plt.figure()
ax = plt.subplot(111)

analytical_solution = Solution()
analytical_solution.color = blue_palette[4]
analytical_solution.ref = True
analytical_solution.name = "Analytical"
analytical_solution.ID = "Analytical"

trnsys_solution = Solution()
trnsys_solution.color = green_palette[4]
trnsys_solution.ref = True
trnsys_solution.name = "TRNSYS"
trnsys_solution.ID = "TRNSYS"

fluent_solution = Solution()
fluent_solution.color = green_palette[3]
fluent_solution.ref = True
fluent_solution.name = "FLUENT"
fluent_solution.ID = "FLUENT"

matlab_solution = Solution()
matlab_solution.color = green_palette[2]
matlab_solution.ref = True
matlab_solution.name = "MATLAB"
matlab_solution.ID = "MATLAB"

ss_solution = Solution()
ss_solution.color = purple_palette[4]
ss_solution.ref = False
ss_solution.name = "Kiva: Steady-State"
ss_solution.ID = "SteadyState"

ade_solution = Solution()
ade_solution.color = red_palette[4]
ade_solution.ref = False
ade_solution.name = "Kiva: ADE"
ade_solution.ID = "ADE"

adi_solution = Solution()
adi_solution.color = red_palette[3]
adi_solution.ref = False
adi_solution.name = "Kiva: ADI"
adi_solution.ID = "ADI"

implicit_solution = Solution()
implicit_solution.color = yellow_palette[4]
implicit_solution.ref = False
implicit_solution.name = "Kiva: Implicit"
implicit_solution.ID = "Implicit"

explicit_solution = Solution()
explicit_solution.color = yellow_palette[3]
explicit_solution.ref = False
explicit_solution.name = "Kiva: Explicit"
explicit_solution.ID = "Explicit"

cn_solution = Solution()
cn_solution.color = yellow_palette[2]
cn_solution.ref = False
cn_solution.name = "Kiva: Crank-Nicolson"
cn_solution.ID = "CrankNicolson"

solutions = [analytical_solution,
             trnsys_solution,
             fluent_solution,
             matlab_solution,
             ss_solution,
             ade_solution,
             adi_solution,
             implicit_solution,
             explicit_solution,
             cn_solution]

cases = ['GC10a',
         'GC30a',
         'GC30b',
         'GC30c',
         'GC60b',
         'GC65b']

values = []
names = []
colors = []

print "Read results..."

for soln in solutions:
    names.append(soln.name)
    colors.append(soln.color)

ind = []

i = 1

for case in cases:
    for soln in solutions:
        values.append(getValue(case,soln))
        ind.append(i)
        i+=1
    i+=1
    
# Bar Chart
print "Create figure..."

width = 1

data = ax.bar(ind, values, width, color=colors)

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

legend = ax.legend(data[:10], names[:10], loc='upper center', ncol=5, bbox_to_anchor=(0.5,-0.05),
                   fancybox=True, shadow=True)


plt.show()
#plt.savefig('../figures/GC10a.pdf')

print "Done."



if __name__ == '__main__':
    pass