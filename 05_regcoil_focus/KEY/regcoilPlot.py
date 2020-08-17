def regcoilPlot(filename):

    # coilsPerHalfPeriod is used only to choose the number of contours to plot for the total current potential:
    coilsPerHalfPeriod = 3

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.io import netcdf
    from scipy.interpolate import interp1d
    import sys, os
    import math

    f = netcdf.netcdf_file(filename,'r',mmap=False)
    nfp = f.variables['nfp'][()]
    ntheta_plasma = f.variables['ntheta_plasma'][()]
    ntheta_coil = f.variables['ntheta_coil'][()]
    nzeta_plasma = f.variables['nzeta_plasma'][()]
    nzeta_coil = f.variables['nzeta_coil'][()]
    nzetal_plasma = f.variables['nzetal_plasma'][()]
    nzetal_coil = f.variables['nzetal_coil'][()]
    theta_plasma = f.variables['theta_plasma'][()]
    theta_coil = f.variables['theta_coil'][()]
    zeta_plasma = f.variables['zeta_plasma'][()]
    zeta_coil = f.variables['zeta_coil'][()]
    zetal_plasma = f.variables['zetal_plasma'][()]
    zetal_coil = f.variables['zetal_coil'][()]
    r_plasma  = f.variables['r_plasma'][()]
    r_coil  = f.variables['r_coil'][()]
    chi2_B = f.variables['chi2_B'][()]
    single_valued_current_potential_thetazeta = f.variables['single_valued_current_potential_thetazeta'][()]
    current_potential = f.variables['current_potential'][()]
    Bnormal_from_plasma_current = f.variables['Bnormal_from_plasma_current'][()]
    Bnormal_from_net_coil_currents = f.variables['Bnormal_from_net_coil_currents'][()]
    Bnormal_total = f.variables['Bnormal_total'][()]
    net_poloidal_current_Amperes = f.variables['net_poloidal_current_Amperes'][()]
    # = f.variables[''][()]

    # lambda used to be called alpha, so try both names for backward-compatibility.
    # Also, we used 'lambdas' instead of 'lambda' to avoid conflict with python's keyword lambda.
    try:
        nlambda = f.variables['nlambda'][()]
        lambdas = f.variables['lambda'][()]
    except:
        nlambda = f.variables['nalpha'][()]
        lambdas = f.variables['alpha'][()]

    # K used to be called J, so try both names for backward-compatibility.
    try:
        chi2_K = f.variables['chi2_K'][()]
        K2 = f.variables['K2'][()]
    except:
        chi2_K = f.variables['chi2_J'][()]
        K2 = f.variables['J2'][()]

    f.close()

    if np.max(np.abs(lambdas)) < 1.0e-200:
        print("lambda array appears to be all 0. Changing it to all 1 to avoid a python error.")
        lambdas += 1

    ########################################################
    # Sort in order of lambda, since for a lambda search (general_option=4 or 5),
    # the output arrays are in the order of the search, which is not so convenient for plotting.
    ########################################################

    permutation = np.argsort(lambdas)
    lambdas = lambdas[permutation]
    chi2_K = chi2_K[permutation]
    chi2_B = chi2_B[permutation]
    Bnormal_total = Bnormal_total[permutation,:,:]
    single_valued_current_potential_thetazeta = single_valued_current_potential_thetazeta[permutation,:,:]
    K2 = K2[permutation,:,:]
    current_potential = current_potential[permutation,:,:]
    if lambdas[-1]>1.0e199:
        lambdas[-1] = np.inf

    ########################################################
    # For 3D plotting, 'close' the arrays in u and v
    ########################################################

    r_plasma  = np.append(r_plasma,  r_plasma[[0],:,:], axis=0)
    r_plasma  = np.append(r_plasma,  r_plasma[:,[0],:], axis=1)
    zetal_plasma = np.append(zetal_plasma,nfp)

    r_coil  = np.append(r_coil,  r_coil[[0],:,:], axis=0)
    r_coil  = np.append(r_coil,  r_coil[:,[0],:], axis=1)
    zetal_coil = np.append(zetal_coil,nfp)

    ########################################################
    # Extract cross-sections of the 3 surfaces at several toroidal angles
    ########################################################

    def getCrossSection(rArray, zetal_old, zeta_new):
        zetal_old = np.concatenate((zetal_old-nfp,zetal_old))
        rArray = np.concatenate((rArray,rArray),axis=0)

        x = rArray[:,:,0]
        y = rArray[:,:,1]
        z = rArray[:,:,2]
        R = np.sqrt(x**2 + y**2)


        ntheta = z.shape[1]
        nzeta_new = len(zeta_new)
        R_slice = np.zeros([nzeta_new,ntheta])
        Z_slice = np.zeros([nzeta_new,ntheta])
        for itheta in range(ntheta):
            interpolator = interp1d(zetal_old, R[:,itheta])
            R_slice[:,itheta] = interpolator(zeta_new)
            interpolator = interp1d(zetal_old, z[:,itheta])
            Z_slice[:,itheta] = interpolator(zeta_new)

        return R_slice, Z_slice

    zeta_slices = np.array([0, 0.25, 0.5, 0.75])*2*np.pi/nfp
    R_slice_plasma, Z_slice_plasma = getCrossSection(r_plasma, zetal_plasma, zeta_slices)
    R_slice_coil, Z_slice_coil = getCrossSection(r_coil, zetal_coil, zeta_slices)

    ########################################################
    # Now make plot of surfaces at given toroidal angle
    ########################################################

    Rmin = R_slice_coil.min()
    Rmax = R_slice_coil.max()
    Zmin = Z_slice_coil.min()
    Zmax = Z_slice_coil.max()

    for whichPlot in range(4):
        plt.figure()
        zeta = zeta_slices[whichPlot]
        plt.plot(R_slice_coil[whichPlot,:], Z_slice_coil[whichPlot,:], 'b.-', label='coil')
        plt.plot(R_slice_plasma[whichPlot,:], Z_slice_plasma[whichPlot,:], 'r.-', label='plasma')
        plt.gca().set_aspect('equal',adjustable='box')
        plt.legend(fontsize='x-small')
        plt.title('zeta='+str(zeta))
        plt.xlabel('R [meters]')
        plt.ylabel('Z [meters]')
        plt.xlim([Rmin,Rmax])
        plt.ylim([Zmin,Zmax])

    ########################################################
    # Pick the lambda values to show in the 2D plots
    ########################################################

#     max_nlambda_for_contour_plots = 16
#     numPlots = min(nlambda,max_nlambda_for_contour_plots)
#     ilambda_to_plot = np.sort(list(set(map(int,np.linspace(1,nlambda,numPlots)))))
#     numPlots = len(ilambda_to_plot)
    max_nlambda_for_contour_plots = 1
    numPlots = 1
    ilambda_to_plot = [-1]

    ########################################################
    # Now make plot of chi^2 over lambda scan
    ########################################################

    plt.figure()
    plt.loglog(chi2_K,chi2_B,'.-r')
    plt.xlabel(r'$\chi^2_K$ [A$^2$]')
    plt.ylabel(r'$\chi^2_B$ [T$^2$ m$^2$]')
    for j in range(numPlots):
        plt.plot(chi2_K[ilambda_to_plot[j]-1], chi2_B[ilambda_to_plot[j]-1],'ob')

    if nlambda>1:
        plt.figure()
        plt.loglog(lambdas,chi2_B,'.-r')
        plt.xlabel('lambda [T$^2$ m$^2$ / A$^2$]')
        plt.ylabel(r'$\chi^2_B$ [T$^2$ m$^2$]')
        for j in range(numPlots):
            plt.plot(lambdas[ilambda_to_plot[j]-1], chi2_B[ilambda_to_plot[j]-1],'ob')

    plt.figure()
    plt.semilogy(lambdas,chi2_B,'.-r')
    plt.xlabel('lambda [T$^2$ m$^2$ / A$^2$]')
    plt.ylabel(r'$\chi^2_B$ [T$^2$ m$^2$]')
    for j in range(numPlots):
        plt.plot(lambdas[ilambda_to_plot[j]-1], chi2_B[ilambda_to_plot[j]-1],'ob')

    if nlambda>1:
        plt.figure()
        plt.loglog(lambdas,chi2_K,'.-r')
        plt.xlabel(r'$\lambda$ [T$^2$ m$^2$ / A$^2$]')
        plt.ylabel(r'$\chi^2_K$ [A$^2$]')
        for j in range(numPlots):
            plt.plot(lambdas[ilambda_to_plot[j]-1], chi2_K[ilambda_to_plot[j]-1],'ob')

    plt.figure()
    plt.semilogy(lambdas,chi2_K,'.-r')
    plt.xlabel('lambda [T$^2$ m$^2$ / A$^2$]')
    plt.ylabel(r'$\chi^2_K$ [A$^2$]')
    for j in range(numPlots):
        plt.plot(lambdas[ilambda_to_plot[j]-1], chi2_K[ilambda_to_plot[j]-1],'ob')

    plt.figure()
    plt.figtext(0.5,0.99,"Blue dots indicate the points in the lambda scan that are plotted in later figures",horizontalalignment='center',verticalalignment='top',fontsize='small')

    ########################################################
    # Prepare for plotting current potential
    ########################################################

    numCols = int(np.ceil(np.sqrt(numPlots)))
    numRows = int(np.ceil(numPlots*1.0/numCols))

    mpl.rc('xtick',labelsize=7)
    mpl.rc('ytick',labelsize=7)

    numContours = 20

    ########################################################
    # Plot single-valued part of current potential
    ########################################################

    for whichPlot in range(numPlots):
        plt.figure()
        plt.contourf(zeta_coil, theta_coil, np.transpose(single_valued_current_potential_thetazeta[ilambda_to_plot[whichPlot]-1,:,:]), numContours)
        plt.colorbar()
        plt.xlabel('zeta',fontsize='x-small')
        plt.ylabel('theta',fontsize='x-small')
        plt.title('lambda='+str(lambdas[ilambda_to_plot[whichPlot]-1]),fontsize='x-small')

    plt.tight_layout()
    plt.figtext(0.5,0.99,"Single-valued part of the current potential",horizontalalignment='center',verticalalignment='top',fontsize='small')

    ########################################################
    # Plot total current potential
    ########################################################

    # Set up contours at appropriate levels
    x = net_poloidal_current_Amperes/nfp
    if abs(x) > np.finfo(float).eps:  # if net_poloidal_current_Ampere presented
        contours = np.linspace(-0.5*x,1.5*x,coilsPerHalfPeriod*2*2,endpoint=False)
    else:
        x = np.max(current_potential) # some new value
        contours = np.linspace(-0.5*x,1.5*x,coilsPerHalfPeriod*2*2,endpoint=False)
    d = contours[1]-contours[0]
    contours = contours + d/2

    zeta_coil_extended = np.concatenate((zeta_coil,[2*np.pi/nfp]))
    theta_coil_extended = np.concatenate((theta_coil,[2*np.pi]))

    for whichPlot in range(numPlots):
        plt.figure()
        data = np.transpose(current_potential[ilambda_to_plot[whichPlot]-1,:,:])
        data_extended = np.concatenate((data,data[0:1,:]),axis=0) # Add the repeated point in theta
        data_extended = np.concatenate((data_extended,data_extended[:,0:1]+x),axis=1) # Add the repeated point in zeta +G
        plt.contourf(zeta_coil_extended, theta_coil_extended, data_extended)
        plt.colorbar()
        plt.xlabel('zeta',fontsize='x-small')
        plt.ylabel('theta',fontsize='x-small')
        plt.title('lambda='+str(lambdas[ilambda_to_plot[whichPlot]-1]),fontsize='x-small')

    plt.figtext(0.5,0.99,"Total current potential",horizontalalignment='center',verticalalignment='top',fontsize='small')

    ########################################################
    # Plot the current density K
    ########################################################

    for whichPlot in range(numPlots):
        plt.figure()
        plt.contourf(zeta_coil, theta_coil, np.sqrt(np.transpose(K2[ilambda_to_plot[whichPlot]-1,:,:])), numContours)
        plt.colorbar()
        plt.xlabel('zeta',fontsize='x-small')
        plt.ylabel('theta',fontsize='x-small')
        plt.title('lambda='+str(lambdas[ilambda_to_plot[whichPlot]-1]),fontsize='x-small')

    plt.figtext(0.5,0.99,"K [A/m]",horizontalalignment='center',verticalalignment='top',fontsize='small')

    ########################################################
    # Plot Bnormal
    ########################################################

    plt.figure()
    plt.contourf(zeta_plasma, theta_plasma, np.transpose(Bnormal_from_plasma_current), numContours)
    plt.colorbar()
    plt.xlabel('zeta',fontsize='x-small')
    plt.ylabel('theta',fontsize='x-small')
    plt.title('From plasma current',fontsize='x-small')

    plt.figure()
    plt.contourf(zeta_plasma, theta_plasma, np.transpose(Bnormal_from_net_coil_currents), numContours)
    plt.colorbar()
    plt.xlabel('zeta',fontsize='x-small')
    plt.ylabel('theta',fontsize='x-small')
    plt.title('From net coil currents',fontsize='x-small')


    for whichPlot in range(numPlots-2):
        plt.figure()
        plt.contourf(zeta_plasma, theta_plasma, np.transpose(Bnormal_total[ilambda_to_plot[whichPlot]-1,:,:]), numContours)
        plt.colorbar()
        plt.xlabel('zeta',fontsize='x-small')
        plt.ylabel('theta',fontsize='x-small')
        plt.title('Total, lambda='+str(lambdas[ilambda_to_plot[whichPlot]-1]),fontsize='x-small')

    plt.figtext(0.5,0.99,"Bnormal [Tesla]",horizontalalignment='center',verticalalignment='top',fontsize='small')

    plt.show()

