import os
import sys
sys.path.append(os.path.abspath('../'))
from riverbuilder.core import river
import matplotlib
import matplotlib.pyplot as plt

#############################################################################################

case_names = ['SFE_Leggett_FB_V']


dir_orig = os.path.abspath('../')

for case_name in case_names:

    dir = os.path.join(dir_orig, 'examples_custom', case_name)
    os.chdir(dir)
    fname = case_name + ".txt"
    outfolder = case_name
    log = case_name + "_log.txt"

    # Run RiverBuilder
    river.buildRiver(fname, outfolder, log)

    plt.close('all')