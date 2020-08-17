def cutCoilsFromRegcoil(filename,nescinFilename,coilsPerHalfPeriod,thetaShift,ilambda):
    import sys

    filename_trimmed = filename.split('/')[-1]
    if filename_trimmed[:12] != 'regcoil_out.':
        print("Error! First argument should be regcoil_out.XXX")
        sys.exit(1)
    coilsFilename='coils.'+filename_trimmed[12:-3]
    print("coilsFilename:",coilsFilename)
    print("coilsPerHalfPeriod:",coilsPerHalfPeriod)
    print("thetaShift:",thetaShift)
    print("ilambda:",ilambda)

    from scipy.io import netcdf
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np

    f = netcdf.netcdf_file(filename,'r',mmap=False)
    theta = f.variables['theta_coil'][()]
    zeta = f.variables['zeta_coil'][()]
    nfp = f.variables['nfp'][()]
    net_poloidal_current_Amperes = f.variables['net_poloidal_current_Amperes'][()]
    current_potential = f.variables['current_potential'][()]
    f.close()

    if abs(net_poloidal_current_Amperes) > np.finfo(float).eps:
        data = current_potential[ilambda,:,:] / net_poloidal_current_Amperes * nfp
    else:
        data = current_potential[ilambda,:,:] / np.max(current_potential[ilambda,:,:])

    # First apply 'roll' to be sure I use the same convention as numpy:
    theta = np.roll(theta,thetaShift)
    # Now just generate a new monotonic array with the correct first value:
    theta = theta[0] + np.linspace(0,2*np.pi,len(theta),endpoint=False)

    data = np.roll(data,thetaShift,axis=1)

    d = 2*np.pi/nfp
    zeta_3 = np.concatenate((zeta-d, zeta, zeta+d))
    data_3 = np.concatenate((data-1,data,data+1))

    fig = plt.figure()
    fig.patch.set_facecolor('white')

    contours = np.linspace(-1,2,coilsPerHalfPeriod*2*3+1)
    d = contours[1]-contours[0]
    contours = contours + d/2
    cdata = plt.contour(zeta_3,theta,np.transpose(data_3),contours)
    plt.colorbar()
    plt.xlabel('zeta')
    plt.ylabel('theta')

    # Repeat with just the contours we care about:
    contours = np.linspace(0,1,coilsPerHalfPeriod*2,endpoint=False)
    d = contours[1]-contours[0]
    contours = contours + d/2
    cdata = plt.contour(zeta_3,theta,np.transpose(data_3),contours,colors='k')

    numCoilsFound = len(cdata.collections)
    if numCoilsFound != 2*coilsPerHalfPeriod:
        print("WARNING!!! The expected number of coils was not the number found.")

    contour_zeta=[]
    contour_theta=[]
    numCoils = 0
    for j in range(numCoilsFound):
        p = cdata.collections[j].get_paths()[0]
        v = p.vertices
        # Make sure the contours have increasing theta:
        if v[1,1]<v[0,1]:
            v = np.flipud(v)

       # close the contours by adding a copy of the first point to the end
        for jfp in range(nfp):
            d = 2*np.pi/nfp*jfp
            contour_zeta.append(v[:,0]+d)
            contour_theta.append(v[:,1])
            numCoils += 1
#             plt.plot(contour_zeta[-1],contour_theta[-1],'.-r',linewidth=1)
#             plt.plot(contour_zeta[-1][0],contour_theta[-1][0],'sk')
          

    # Now read nescin filename to map the theta-zeta coordinates of the contours to xyz coordinates.
    f=open(nescinFilename,'r')

    line = ''
    while "np     iota_edge       phip_edge       curpol" not in line:
        line = f.readline()
    line = f.readline()
    nfp_nescin = int(line.split()[0])
    if nfp != nfp_nescin:
        print("Error! nfp from regcoil_out does not match nfp from nescin!")
        sys.exit(1)

    contour_R = []
    contour_Z = []
    for j in range(numCoils):
        contour_R.append(contour_zeta[j]*0)
        contour_Z.append(contour_zeta[j]*0)

    line = ''
    #while "Number of fourier modes in table" not in line:
    while "------ Current Surface" not in line:
        line = f.readline()
    line = f.readline()
    line = f.readline()
    print("Number of Fourier modes in coil surface from nescin file: ",line)
    nmodes = int(line)
    line = f.readline()
    line = f.readline()
    for imode in range(nmodes):
        data = f.readline().split()
        m = int(data[0])
        n = -int(data[1]) * nfp
        # Sign flip in n because bnormal uses NESCOIL convention.
        # See bn_read_vmecf90.f line 89.
        crc = float(data[2])
        czs = float(data[3])
        crs = float(data[4])
        czc = float(data[5])
        # Skip remaining columns
        for j in range(numCoils):
            angle = m*contour_theta[j] - n*contour_zeta[j]
            contour_R[j] = contour_R[j] + crc*np.cos(angle) + crs*np.sin(angle)
            contour_Z[j] = contour_Z[j] + czs*np.sin(angle) + czc*np.cos(angle)
    f.close()

    contour_X = []
    contour_Y = []
    maxR=0
    for j in range(numCoils):
        maxR = np.max((maxR,np.max(contour_R[j])))
        contour_X.append(contour_R[j]*np.cos(contour_zeta[j]))
        contour_Y.append(contour_R[j]*np.sin(contour_zeta[j]))

    coilCurrent = net_poloidal_current_Amperes / numCoils

    # Find the point of minimum separation
    minSeparation2=1.0e+20
    for whichCoil1 in range(numCoils):
        for whichCoil2 in range(whichCoil1):
            for whichPoint in range(len(contour_X[whichCoil1])):
                dx = contour_X[whichCoil1][whichPoint] - contour_X[whichCoil2]
                dy = contour_Y[whichCoil1][whichPoint] - contour_Y[whichCoil2]
                dz = contour_Z[whichCoil1][whichPoint] - contour_Z[whichCoil2]
                separation2 = dx*dx+dy*dy+dz*dz
                this_minSeparation2 = np.min(separation2)
                if this_minSeparation2<minSeparation2:
                    minSeparation2 = this_minSeparation2
                    x1 = contour_X[whichCoil1][whichPoint]
                    y1 = contour_Y[whichCoil1][whichPoint]
                    z1 = contour_Z[whichCoil1][whichPoint]
                    index=np.argmin(separation2)
                    x2 = contour_X[whichCoil2][index]
                    y2 = contour_Y[whichCoil2][index]
                    z2 = contour_Z[whichCoil2][index]
                
    print("Minimum coil separation:",np.sqrt(minSeparation2))

    # Write coils file
    f = open(coilsFilename,'w')
    f.write('periods '+str(nfp)+'\n')
    f.write('begin filament\n')
    f.write('mirror NIL\n')

    for j in range(numCoils):
        N = len(contour_X[j])
        for k in range(N):
            f.write('{:14.22e} {:14.22e} {:14.22e} {:14.22e}\n'.format(contour_X[j][k],contour_Y[j][k],contour_Z[j][k],coilCurrent))
        # Close the loop
        k=0
        f.write('{:14.22e} {:14.22e} {:14.22e} {:14.22e} 1 Modular\n'.format(contour_X[j][k],contour_Y[j][k],contour_Z[j][k],0))

    f.write('end\n')
    f.close()

    plt.show()
