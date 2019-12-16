import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

## data directories need to be of the form
## Parallel test type
## -- all four different threaded files for the test type
## Parallel test type
## -- all four different threaded files for the test type

## I do not know how it orders the files on linux. But on windows I follow the order that it iterates through
## the files when dictating what line corresponds to what title in the legend



legend_list = ["One Thread", "Two Threads", "Three Threads", "Four Threads"]
lines = ['solid', "dashed", "dashdot", "dotted"]
types = ["Jacobi", "Gauss_Seidel"]
cur_dir = os.getcwd()
dirs = [os.path.join(cur_dir, "jacobi_OMP"), os.path.join(cur_dir, "testGauss_OMP")]
fig = plt.figure()
list_of_data = []
dim = np.arange(10,301,1)
test_type_counter = 0

for dir in dirs:

    fig = plt.figure()
    list_of_data = []
    line_type_counter = 0
    
    for f in os.listdir(dir):
        print(f)
        l = []
        f = open(os.path.join(dir, f), "r")
        if f.mode == "r":
            content = f.read()
            
            content = content.split()
        
            for item in content:
                l.append(float(item)/1000000)
                
                list_of_data.append(l.copy())
                
                fig.suptitle('Speed Up Between Thread Counts ' + types[test_type_counter] )
                plt.xlabel('N', fontsize=18)
                plt.ylabel('Time (Seconds)', fontsize=16)
                
                for i in range(len(list_of_data)):
                    plt.plot(dim, list_of_data[i], linestyle=lines[line_type_counter])
                    line_type_counter += 1
                    line_type_counter = line_type_counter % 4

                plt.legend(legend_list)
                fig.savefig(types[test_type_counter] + ".jpg")
                
                plt.close(fig)
                list_of_data = []
                time_fig = plt.figure()
                line_type_counter = 0
                test_type_counter += 1
                
                #smoothed = gaussian_filter1d(l, sigma=2)
                #plt.plot(smoothed)
                                                                                                                                                                                                
