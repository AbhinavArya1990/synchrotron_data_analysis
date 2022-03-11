# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 17:16:41 2017

"""
#==============================================================================
# #import functions
#==============================================================================
import warnings
warnings.filterwarnings("ignore")
import pyFAI
import fabio
import numpy
import os
import glob
import matplotlib.pyplot as plt
#==============================================================================
# # Input
#==============================================================================
manual = 0
find_azt = 1
save_img_set = 0
plot_int= 0
display = 0
loop = 0      # 0 for NONE (1 for RANGE, 2 for two set of angles)
twoD = 0
oneD = 0
save = 0
# get all files
file_path = "C:\Users\puruschr\Desktop\First beamtime Synchrotron PSI December 2016\pyFAI-test\calibration\EP14\\"
all_files = glob.glob(file_path+"*.cbf")
all_files.sort()
print "number of files found: "+ str(len(all_files))
# radial and azimuth range
if loop == 0 and save_img_set == 0:
    print 'no range of azimuth inputed'
else:
    az11 = float(raw_input('Enter starting azimuth range in degree (smallest/ negative): '))
    az21 = float(raw_input('Enter ending azimuth range in degree (largest): '))
    step = int(raw_input('Enter number of angles needed (integer): '))
    az_range = numpy.linspace(az11, az21, num=step)
#==============================================================================
# # creating multi integration scheme (INPUT from PONI- Manual)
#==============================================================================
if manual == 1:
    Distance= 1.01778467929
    PixelSize1= 0.000172
    PixelSize2= 0.000172
    Poni1= 0.249736153863
    Poni2= 0.461275616348
    Rot1= -0.173505121333
    Rot2= 0.00264538851228
    Rot3= -3.141592653589793
    wl = 4.76977855865e-11
    # getting detector from dictionary
    dete = pyFAI.detectors.Pilatus6M()
    print(dete)
    # generating geometry
    ai = pyFAI.AzimuthalIntegrator(dist=Distance, poni1=Poni1,\
             poni2=Poni2, detector=dete, rot1 = Rot1, rot2= Rot2, rot3 =Rot3, wavelength = wl)
    print(ai)
#==============================================================================
# # Automatic input from PONI file
#==============================================================================
elif manual == 0:
    ai = pyFAI.load("C:\Users\puruschr\Desktop\First beamtime Synchrotron PSI December 2016\pyFAI-test\calibration\LaB6_26keV_ff_0002.poni")
    print(ai)
#==============================================================================
# # Find azimuthal range
#==============================================================================
if find_azt == 1:
    img_file = all_files[25]
    img1 = fabio.open(img_file).data 
    # 2D integration
    az1, az2 = -40, 40
    I1, tth1, chi1 = ai.integrate2d(img1, 300, 360, azimuth_range=(az1, az2), unit="2th_deg",\
                                            method="BBox") #azimuthal integration
    """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr","""
    tth, I = ai.integrate1d(img1, 2464, unit="2th_deg", azimuth_range=(az1, az2),\
                                            )
    """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr", "nosplit_csr", "full_csr","""

    ## PLOT
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(2,2,1)
    ax.imshow(numpy.log(img1), origin="lower")
    ax1 = plt.subplot(2,2,2)
    ax1.imshow(numpy.log(I1), origin="lower", extent=[tth1.min(), tth1.max(), chi1.min(), chi1.max()], aspect="auto")
    ax1.set_xlabel("2 theta (deg)")
    ax1.set_ylabel("Azimuthal angle chi (deg)")
    ax2 = plt.subplot(2,2,3)
    ax2.imshow(numpy.log(img1), origin="lower")
    ax3 = plt.subplot(2,2,4,sharex=ax1)
    ax3.scatter(tth, I)
    ax3.set_xlabel("2 theta (deg)")
    plt.show()

if save_img_set == 1:
    for i in range(len(az_range)-1):
        az1 = az_range[i]
        az2 = az_range[i+1]
        img_file = all_files[0]
        img1 = fabio.open(img_file).data 
        newFilename2Dimg = file_path+"2D_"+str(az1)+"_"+str(az2)+"_.jpg"
        # 2D integration
        I1, tth1, chi1 = ai.integrate2d(img1, 300, 360, azimuth_range=(az1, az2), unit="2th_deg",\
                                                method="BBox") #azimuthal integration
        """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr","""
        tth, I = ai.integrate1d(img1, 2464, unit="2th_deg", azimuth_range=(az1, az2),\
                                                method = "BBox")
        """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr", "nosplit_csr", "full_csr","""
        
        ## PLOT
        fig = plt.figure(figsize=(10,10))
        ax = plt.subplot(2,2,1)
        ax.imshow(numpy.log(img1), origin="lower")
        ax1 = plt.subplot(2,2,2)
        ax1.imshow(numpy.log(I1), origin="lower", extent=[tth1.min(), tth1.max(), chi1.min(), chi1.max()], aspect="auto")
        ax1.set_xlabel("2 theta (deg)")
        ax1.set_ylabel("Azimuthal angle chi (deg)")
        ax2 = plt.subplot(2,2,3)
        ax2.imshow(numpy.log(img1), origin="lower")
        ax3 = plt.subplot(2,2,4,sharex=ax1)
        ax3.scatter(tth, I)
        ax3.set_xlabel("2 theta (deg)")
        fig.savefig(newFilename2Dimg, bbox_inches='tight',format='jpg', dpi=1000) 
        plt.close(fig)
        print 'Azimuth range '+str(az1)+ ' to '+str(az2)+ ' completed'

if plot_int == 1:
    img_file = all_files[10]
    img1 = fabio.open(img_file).data 
    # 2D integration
    xplot, xplot1= [], []
    yplot, yplot1= [], []
    az_range = numpy.linspace(-40, 40, num=160)
    for i in range(len(az_range)-1):
        az1 = az_range[i]
        az2 = az_range[i+1]
        I1, tth1, chi1 = ai.integrate2d(img1, 300, 360, azimuth_range=(az1, az2), unit="2th_deg",\
                                                method="BBox" , dummy=-1) #azimuthal integration
        """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr","""
        tth, I = ai.integrate1d(img1, 2464, unit="2th_deg", azimuth_range=(az1, az2),\
                                                 dummy=-1)
        """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr", "nosplit_csr", "full_csr","""
        idx = (tth>13)*(tth<14) # for 111 reflection
        ind = numpy.where(idx)
        intensity = I[ind]
        xplot.append(numpy.linspace(az1, az1, num=len(ind[0])))
        yplot.append(intensity)
        xplot1.append(az1)
        yplot1.append(numpy.max(intensity))
        print 'Azimuth range '+str(az1)+ ' to '+str(az2)+ ' completed'
    ## PLOT
    fig = plt.figure(figsize=(10,10))
    plt.scatter(xplot1, yplot1)
    find_max = numpy.where(yplot1==numpy.max(yplot1))
    plt.axvline(x=xplot1[find_max[0][0]],linestyle='dashed')
    plt.annotate('maximum intensity at azimuthal angle of '+str(xplot1[find_max[0][0]]), \
                xy=(xplot1[find_max[0][0]], numpy.max(yplot1)), \
                xytext=(xplot1[find_max[0][0]]+1, numpy.max(yplot1)+1),)
    plt.xlabel("Azimuthal angle chi (deg)")
    plt.ylabel("Intensity")
    fig.savefig('azi_range.jpg', bbox_inches='tight',format='jpg', dpi=1000) 
    plt.close(fig)     
    
    fig = plt.figure(figsize=(10,10))
    plt.scatter(xplot, yplot)
    plt.xlabel("Azimuthal angle chi (deg)")
    plt.ylabel("Intensity")
    fig.savefig('azi_range1.jpg', bbox_inches='tight',format='jpg', dpi=1000) 
    plt.close(fig)
#==============================================================================
# multi file analysis
#==============================================================================
file_path = file_path + "dat_files\\"
if not os.path.exists(file_path):
    os.makedirs(file_path)
    # create a directory for saving results
if loop == 1:
    for i in range(len(az_range)-1):
        az1 = az_range[i]
        az2 = az_range[i+1]
        file_no = 1 #initial file number
        for one_file in all_files:
            destination = os.path.splitext(one_file)[0]+".dat"
            image = fabio.open(one_file).data
            newFilename2D = file_path+"2D_az_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".dat"
            newFilename1D = file_path+"1D_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".dat"
            newFilename2Dimg = file_path+"2D_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".jpg"
            newFilename1Dimg = file_path+"1D_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".jpg"
            
            if twoD == 1:
                # Do 2D integration
        #        I, tth, chi = ai.integrate2d(image, 300, 360, azimuth_range=(140,180), \
        #                            unit="2th_deg", dummy=-1, delta_dummy=1, filename=newFilename2D) #azimuthal integration
                I, tth, chi = ai.integrate2d(image, 300, 360, azimuth_range=(az1, az2), \
                                    unit="2th_deg", method = "BBox", filename=newFilename2D)        
                """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr","""
                if display == 1:            
                    ## PLOT
                    fig = plt.figure(figsize=(10,10))
                    ax = plt.subplot(1,2,1)
                    ax.imshow(numpy.log(image), origin="lower")
                    ax = plt.subplot(1,2,2)
                    ax.imshow(I, origin="lower", extent=[tth.min(), tth.max(), chi.min(), chi.max()], aspect="auto")
                    ax.set_xlabel("2 theta (deg)")
                    ax.set_ylabel("Azimuthal angle chi (deg)")  
                    plt.show()
                    
                if save == 1:
                    fig = plt.figure(figsize=(10,10))
                    ax = plt.subplot(1,2,1)
                    ax.imshow(numpy.log(image), origin="lower")
                    ax = plt.subplot(1,2,2)
                    ax.imshow(I, origin="lower", extent=[tth.min(), tth.max(), chi.min(), chi.max()], aspect="auto")
                    ax.set_xlabel("2 theta (deg)")
                    ax.set_ylabel("Azimuthal angle chi (deg)") 
                    fig.savefig(newFilename2Dimg, bbox_inches='tight',format='jpg', dpi=1000) 
                    plt.close(fig)
            # Do 1D integration
            if oneD == 1:  
        #        tth, I = ai.integrate1d(image, 2500, unit="2th_deg", azimuth_range=(140,180),\
        #                                            filename=newFilename1D, dummy=-1, delta_dummy=1)
                tth, I = ai.integrate1d(image, 2464, unit="2th_deg", azimuth_range=(az1, az2),\
                                                    filename=newFilename1D)#, method = "csr")
                """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr", "nosplit_csr", "full_csr","""
                if display == 1:            
                    ## PLOT
                    fig = plt.figure(figsize=(10,10))
                    ax = plt.subplot(1,2,1)
                    ax.imshow(numpy.log(image), origin="lower")
                    ax = plt.subplot(1,2,2)
                    ax.plot(tth, I)
                    ax.set_xlabel("2 theta (deg)")
                    plt.show()
                
                if save == 1:
                    fig = plt.figure(figsize=(10,10))
                    ax = plt.subplot(1,2,1)
                    ax.imshow(numpy.log(image), origin="lower")
                    ax = plt.subplot(1,2,2)
                    ax.plot(tth, I)
                    ax.set_xlabel("2 theta (deg)")
                    fig.savefig(newFilename1Dimg, bbox_inches='tight',format='jpg', dpi=1000) 
                    plt.close(fig)
            
            file_no += 1 
        print 'Azimuth range '+str(az1)+ ' to '+str(az2)+ ' completed'

elif loop == 2:
#    az1 = az11
#    az2 = az21
    az1 = -0.2    # smallest number
    az2 = 3.3
    file_no = 1 #initial file number
    for one_file in all_files:
        destination = os.path.splitext(one_file)[0]+".dat"
        image = fabio.open(one_file).data
        newFilename2D = file_path+"2D_az_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".dat"
        newFilename1D = file_path+"1D_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".dat"
        newFilename2Dimg = file_path+"2D_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".jpg"
        newFilename1Dimg = file_path+"1D_"+str(az1)+"_"+str(az2)+"_"+str(file_no)+".jpg"
        
        if twoD == 1:
            # Do 2D integration
    #        I, tth, chi = ai.integrate2d(image, 300, 360, azimuth_range=(140,180), \
    #                            unit="2th_deg", dummy=-1, delta_dummy=1, filename=newFilename2D) #azimuthal integration
            I, tth, chi = ai.integrate2d(image, 300, 360, azimuth_range=(az1, az2), \
                                unit="2th_deg", method = "BBox", dummy=-1, filename=newFilename2D)        
            """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr","""
            if display == 1:            
                ## PLOT
                fig = plt.figure(figsize=(10,10))
                ax = plt.subplot(1,2,1)
                ax.imshow(numpy.log(image), origin="lower")
                ax = plt.subplot(1,2,2)
                ax.imshow(I, origin="lower", extent=[tth.min(), tth.max(), chi.min(), chi.max()], aspect="auto")
                ax.set_xlabel("2 theta (deg)")
                ax.set_ylabel("Azimuthal angle chi (deg)")  
                plt.show()
                
            if save == 1:
                fig = plt.figure(figsize=(10,10))
                ax = plt.subplot(1,2,1)
                ax.imshow(numpy.log(image), origin="lower")
                ax = plt.subplot(1,2,2)
                ax.imshow(I, origin="lower", extent=[tth.min(), tth.max(), chi.min(), chi.max()], aspect="auto")
                ax.set_xlabel("2 theta (deg)")
                ax.set_ylabel("Azimuthal angle chi (deg)") 
                fig.savefig(newFilename2Dimg, bbox_inches='tight',format='jpg', dpi=1000) 
                plt.close(fig)
        # Do 1D integration
        if oneD == 1:  
    #        tth, I = ai.integrate1d(image, 2500, unit="2th_deg", azimuth_range=(140,180),\
    #                                            filename=newFilename1D, dummy=-1, delta_dummy=1)
            tth, I = ai.integrate1d(image, 2464, unit="2th_deg", azimuth_range=(az1, az2),\
                                                filename=newFilename1D, dummy=-1)
            """method: can be "numpy", "cython", "BBox" or "splitpixel", "lut", "csr", "nosplit_csr", "full_csr","""
            if display == 1:            
                ## PLOT
                fig = plt.figure(figsize=(10,10))
                ax = plt.subplot(1,2,1)
                ax.imshow(numpy.log(image), origin="lower")
                ax = plt.subplot(1,2,2)
                ax.plot(tth, I)
                ax.set_xlabel("2 theta (deg)")
                plt.show()
            
            if save == 1:
                fig = plt.figure(figsize=(10,10))
                ax = plt.subplot(1,2,1)
                ax.imshow(numpy.log(image), origin="lower")
                ax = plt.subplot(1,2,2)
                ax.plot(tth, I)
                ax.set_xlabel("2 theta (deg)")
                fig.savefig(newFilename1Dimg, bbox_inches='tight',format='jpg', dpi=1000) 
                plt.close(fig)
        
        file_no += 1 
    print 'Azimuth range '+str(az1)+ ' to '+str(az2)+ ' completed'