print "Import libraries..."

from gen_bestest import *

def getValues(solution):
    if (solution.ref):
        return readXLS(solution.ID)
    else:
        return readTimeseries(solution.ID)

def readXLS(ID):
    worksheet = workbook.sheet_by_name(ID)
    if (ID == 'FLUENT'):
        col = 5
    else:
        col = 4
    
    return worksheet.col_values(col,130,8890)

def readTimeseries(ID):
    df = pd.read_csv('../GC40a/' + ID + '-IA/Timeseries.csv',
                    header=0,
                    names=['time','W'],
                    parse_dates=True,
                    index_col=0)
    
    return list(df.ix['2103'].W)
        
print "Read results..."    

solutions = []
solutions.append(trnsys_solution)
solutions.append(fluent_solution)
solutions.append(matlab_solution)
#solutions.append(ade_solution)
solutions.append(implicit_solution)

# Line Chart
print "Create figure..."
file_name = 'bestest_gc40a'

# Figure style
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_context("paper", {'axes.labelsize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

fig = plt.figure()
ax = fig.add_subplot(111)

for soln in solutions:
    print "...Plotting: " + soln.name
    soln.values = getValues(soln)
    ax.plot(soln.values, label=soln.name)


ax.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5,-0.05), fancybox=True)
ax.yaxis.grid()

plt.show()
#fig.savefig(output_dir + 'images/' + file_name + '.pdf')

print "Done."



if __name__ == '__main__':
    pass