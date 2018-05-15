#!/usr/bin/python
import average_domain_size as ave_dom_sz
######################################################################
# This script takes three commandline arguments:
#
#   argument1: the fullpath name of the directory with the data files
#   argument2: the name of the directory in which to save to ads plot
#   argument3: a string to be used as the plot filename
#
######################################################################

############################################################################
# main script
############################################################################

data_directory = argv[1]
save_directory = argv[2]
plot_name = argv[3]

filelist = os.listdir(data_directory)
filelist = fnmatch.filter(filelist,'c1*')
filelist.sort(key=lambda x: int(x.split('.')[0].split('_')[1]))

# create ads array
points=len(filelist)
ads = np.zeros(points)
i = 0

# calculate an average domain size (ads) for each data file
for filename in filelist:
    ads[i] = ave_dom_sz.calc_ads(data_directory+'/'+filename)
    i = i + 1

# save data
skip = int(filelist[points-1].split('.')[0].split('_')[1])/(points-1)
t = np.arange(0,points*skip,skip) 

# save the ads data for replotting or manipulation later
np.savez(save_directory+'/'+plot_name+'.npz',ads,t)
