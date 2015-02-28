######################################################################
##* Script to plot data
##* June 2014
##* Sarah Shipley
######################################################################

import pandas as pd
import matplotlib.pyplot as plt

compound = raw_input('ROOT: ')
phase = raw_input('gas or aer phase? ')
if phase == 'gas':
    key = compound
elif phase == 'aer':
    key = compound + '_' + 'aer'
else:
    print 'Invalid Response!'
    key = 'none'

dframe = pd.DataFrame.from_csv(key+'_formatted.txt')

#plt.plot(dframe.index,dframe.APINENE)
plt.plot(dframe.index,dframe)
plt.ion
plt.show()


