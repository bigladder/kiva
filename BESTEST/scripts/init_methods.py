print "Import libraries..."

from gen_bestest import *

years = ['2097','2098','2099','2100','2101','2102','2103']
        
class Result:
    def __init__(self, methodID):
        self.methodID = methodID
        if methodID == 'C':
            self.method = "Constant"
        elif methodID == 'K':
            self.method = "Kusuda"
        elif methodID == 'SS':
            self.method = 'Steady-State'
        else: # methodID == 'IA':
            self.method = 'Implicit Acceleration'
        
        self.df = pd.read_csv('../GC40a/Implicit-' + methodID + '/Timeseries.csv',
                              header=0,
                              names=['time','W'],
                              parse_dates=True,
                              index_col=0)

        sums = []
        diffs = []
        for year in years:
            sum = self.df.ix[year].W.sum()
            sums.append(sum)
            
        for i in range(0,7):
            if i == 0:
                diffs.append(None)
            else:
                diff = abs(sums[i] - sums[i-1])/sums[i-1]
                diffs.append(diff)
            
        self.sums = sums
        self.diffs = diffs
        
        

print "Read results..."    

results = {}
methods = ['C','SS','IA']

for method in methods:
    print "...Reading: " + method
    results[method] = Result(method)


# Line Chart
print "Create figure..."

file_name = "init_methods"

# Figure style
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_context("paper", {'axes.labelsize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

year_1_fig = plt.figure()
year_1_ax = year_1_fig.add_subplot(111)


for method in methods:
    print "...Plotting: " + method
    #year_1_ax.plot(list(results[method].df.ix['2097'].W), label=results[method].method)
    year_1_ax.plot(list(results[method].df.W), label=results[method].method)

#print "...Plotting: Final"
#year_1_ax.plot(list(results['IA'].df.ix['2103'].W), label='Final')

box = year_1_ax.get_position()
year_1_ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5,-0.05), fancybox=True)
year_1_ax.yaxis.grid()
year_1_ax.set_ylim([1500,5000])

plt.show()
#year_1_fig.savefig(output_dir + 'images/' + file_name + '.pdf')

for method in methods:
    print "...Diffs: " + method
    print results[method].diffs


print "Done."



if __name__ == '__main__':
    pass