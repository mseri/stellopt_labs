
def mangoPlot(mango_filenames):
# Copyright 2019, University of Maryland and the MANGO development team.
#
# This file is part of MANGO.
#
# MANGO is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# MANGO is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MANGO.  If not, see
# <https://www.gnu.org/licenses/>.

    myfigsize=(14,6.8)

    marker_size = 2
    line_width = 0.7

    import os
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    import glob

    files_function_evaluations = []
    files_objective_function = []
    files_times = []
    filenames = []
    files_parameters = []

    for j in range(len(mango_filenames)):
        filename = mango_filenames[j]
        if os.path.isfile(filename):
            filenames.append(filename)
            basename = os.path.basename(filename)
            if basename[:9] != "mango_out":
                print("WARNING: Including file "+filename+" even though it does not begin with mango_out")

    print("Files that will be read and plotted:")
    for file in filenames:
        print("  "+file)
    print()
    
    for k in range(len(filenames)):
        filename = filenames[k]
        f = open(filename,'r')
        lines = f.readlines()
        f.close()

        temp = lines[3].split(',')
        try:
            N_parameters = int(temp[0])
        except:
            print("ERROR! Unable to read N_parameters from line 3 of "+filename)
            print("This probably means this file is not a correctly formatted mango_out file.")
            raise

        function_evaluations = []
        times = []
        objective_function = []
        parameters = np.zeros((1,N_parameters))
        for j in range(5,len(lines)):
            temp = lines[j].split(',')
            try:
                function_evaluations.append(int(temp[0]))
            except:
                print("ERROR! Unable to convert "+temp[0]+" to int on line "+str(j)+" of file "+filename)
                print("This probably means this file is not a correctly formatted mango_out file.")
                raise

            try:
                times.append(float(temp[1]))
            except:
                print("ERROR! Unable to convert "+temp[1]+" to float on line "+str(j)+" of file "+filename)
                print("This probably means this file is not a correctly formatted mango_out file.")
                raise
                
            try:
                if (j == 5):
                    parameters[0,:] = temp[2:N_parameters+2]
                else:
                    parameters = np.vstack((parameters,temp[2:N_parameters+2]))
            except:
                print("ERROR! Unable to convert "+str(temp[2:N_parameters+2])+" to float on line "+str(j)+" of file "+filename)
                print("This probably means this file is not a correctly formatted mango_out file.")
                raise
                
            try:
                this_objective_function = float(temp[N_parameters+2])
            except:
                print("Warning: unable to convert "+temp[N_parameters+2]+" to float in file "+filename)
                this_objective_function = np.nan

            # Stellopt sets failed results to 1e+12, which makes it hard to see the interesting structure in the objective function for successful runs.
            # So let's just not show failed runs.
            if this_objective_function > 1.0e+11:
                this_objective_function = np.nan


            objective_function.append(this_objective_function)

        if k==0:
            min_objective_function = np.nanmin(objective_function)
            # np.nanmin is a minimum excluding any nans.
        else:
            if len(objective_function) > 0: # Failed runs have len(objective_function)=0, causing np.nanmin to fail
                min_objective_function = np.nanmin((min_objective_function, np.nanmin(objective_function)))

        files_function_evaluations.append(function_evaluations[:-1])
        files_times.append(times[:-1])
        files_objective_function.append(objective_function[:-1])
        files_parameters.append(parameters[:-1])

    N_files = len(files_function_evaluations)
    print("Minimum objective function found:",min_objective_function)

    #########################################################
    # Done reading files. Now make the plot.
    #########################################################

    fig = plt.figure(figsize=myfigsize)
    fig.patch.set_facecolor('white')

    numCols = 1
    numRows = 3
    plotNum = 1

    linespecs = ['o','^','s','v']

    plt.figure()
    for j in range(N_files):
        linespec=linespecs[np.mod(int(np.floor(j/10)),len(linespecs))]+'-'
        plt.semilogy(files_function_evaluations[j], files_objective_function[j] - min_objective_function, linespec, label = filenames[j], markersize=marker_size, linewidth=line_width)
    plt.xlabel('Function evaluation')
    plt.ylabel('(Objective function) - (min objective function)')
    plt.grid(True)
    ncol=int(np.floor(N_files/20))+1
    if N_files > 1:
        plt.legend(loc='upper right',fontsize=7,ncol=ncol)
      
    for k in range(len(files_parameters[0][0,:])):
        plt.figure()
#         ax = plt.gca()
#         ax.set_yticks(ax.get_yticks()[::2])
        plt.title('Parameter '+str(k+1))
        for j in range(N_files):
#             linespec=linespecs[np.mod(int(np.floor(j/10)),len(linespecs))]+'-'
            plt.plot(files_function_evaluations[j], files_parameters[j][:,k],label = filenames[j], markersize=marker_size, linewidth=line_width)
            plt.xlabel('Function evaluation')
            plt.ylabel('Parameter value')
            ncol=int(np.floor(N_files/20))+1
            if N_files > 1:
                plt.legend(loc='upper right',fontsize=7,ncol=ncol) 
            ax = plt.gca()
            ax.yaxis.set_major_locator(plt.MaxNLocator(10))

    ##############################################################

    plt.show()

