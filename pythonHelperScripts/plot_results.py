import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

## data directories need to be of the form
## test type
## ----1
## --------iterations
## --------time
## ----2
## --------iterations
## --------time
## ----3
## --------iterations
## --------time
## ----4
## --------iterations
## --------time
## test type
## ....

## I do not know how it orders the files on linux. But on windows I follow the order that it iterates through
## the files when dictating what line corresponds to what title in the legend



legend_list = ["Gauss - Iterative", "Gauss - Parallel Ring", "Gauss - Ring", "Jacobi - Ring", "Jacobi - Parallel Ring", "Jacobi - Iterative"]
lines = ['solid', "dashed", "dashdot"]
cur_dir = os.getcwd()
dirs = [os.path.join(cur_dir, "1"), os.path.join(cur_dir, "2"), os.path.join(cur_dir, "3"), os.path.join(cur_dir, "4")]
fig = plt.figure()
list_of_data = []
dim = np.arange(10,101,1)


for dir in dirs:

    it_fig = plt.figure()
    iteration_dir = os.path.join(dir, "iterations")
    time_dir = os.path.join(dir, "time")
    list_of_data = []
    j = 0
    
    for f in os.listdir(iteration_dir):
        print(f)
        l = []
        f = open(os.path.join(iteration_dir, f), "r")
        if f.mode == "r":
            content = f.read()
            content = content.split()
        
            for item in content:
                l.append(int(item))
            
                list_of_data.append(l.copy())
                
                it_fig.suptitle('Iteration Comparisons on Relaxation Mesh of Size NxN')
                plt.xlabel('N', fontsize=18)
                plt.ylabel('Iterations', fontsize=16)
                
            for i in range(len(list_of_data)):
                plt.plot(dim, list_of_data[i], linestyle=lines[j])
                j += 1
                j = j % 3

            plt.legend(legend_list)
            it_fig.savefig(dir + "_iterations.jpg")
            
            plt.close(it_fig)
            list_of_data = []
            time_fig = plt.figure()
            j = 0
            
            for f in os.listdir(time_dir):
                print(f)
                l = []
                f = open(os.path.join(time_dir, f), "r")
                if f.mode == "r":
                    content = f.read()
                    content = content.split()
                    
                    for item in content:
                        l.append((float(item)/1000000))

                    list_of_data.append(l.copy())
                    plt.legend(legend_list)
                    time_fig.suptitle('Time Comparisons on Relaxation Mesh of Size NxN')
                    plt.xlabel('N', fontsize=18)
                    plt.ylabel('Time (Seconds)', fontsize=16)
                    for i in range(len(list_of_data)):
                        plt.plot(dim, list_of_data[i], linestyle = lines[j])
                        j += 1
                        j = j % 3

                    plt.legend(legend_list)
                    time_fig.savefig(dir + "_time.jpg")
                    
                    plt.close(fig)

                    #smoothed = gaussian_filter1d(l, sigma=2)
                    #plt.plot(smoothed)
                                                                                                                                                                                                                                                                                                                                                                
