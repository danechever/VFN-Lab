import sys
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.pylab as pl
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm, Normalize
import matplotlib as mpl
from astropy.io import fits, ascii
from astropy.modeling import models
from astropy import modeling
from astropy import units as u
import os
from scipy import signal
import scipy.optimize as opt
import scipy
import multiprocessing as mp
# from FEU_Fiber_cmds import FEU_Fiber_cmds
from Star_Tracker_cmds import Tracking_cmds
from Track_Cam_cmds import TC_cmds
from FAM_cmds import FAM_cmds
# from FIU_Commands import get_path

# from PD_cmds import redPM_cmds
# from FEU_FAM_cmds import FEU_FAM_cmds
# import Acquisition
# import ktl
# import Organization
from PSF_finder import moments, Gaussian2D, PSF_Finder
import time
import warnings
import datetime
# import pdb
from scipy.signal import medfilt as medfilt
sys.path.append('/home/nfiudev/dev/')
from redPM_cmds import redPM_cmds
sys.path.append('/home/nfiudev/dev/dechever/DMCharac/Code/')
from DM import DM

cam = TC_cmds()
track = Tracking_cmds()
FAM = FAM_cmds()
DM = DM()
pd_sensitivity = 10
nreads = 10

nolls = {4:'Defocus', 5:'Astig 45', 6:'Astig 0', 7:'Coma Y', 8:'Coma X', 9:'Trefoil Y', 10:'Trefoil X', 11:'Spherical',
         12:'', 13:'', 14:'', 15:'', 16:'', 17:'', 18:'', 19:'', 20:'', 21:'', 22:'', 23:'', 24:''}

def _get_cred2_flux(plot=False, nb_im=10):
    '''Returns average total flux in a 100x100 cred2 pixel square around the current tracking goal

    Keyword arguments:
    plot: (boolean, default false) Whether to show the cutout
    nb_im: (int, default 10) Number of cred2 frames to average

    Returns:
    Average flux (total counts/integration time) of the cutouts.
    '''
    A = np.zeros(nb_im)
    for i in range(A.size):
        ### Get the processed cred2 image
        procim = track.Proc_im.get_data(reform=True)
        ### Get location of psf
        goal = track.get_psf_goal()
        ### Check SHM tha tthese are right
        xcen = int(np.round(goal[1]))
        ycen = int(np.round(goal[2]))

        ### cutout - size tbd
        cutim = procim[ycen-50:ycen+51, xcen-50:xcen+51]
        if plot:
            plt.imshow(cutim, extent=[xcen-50, xcen+50, ycen-50, ycen+50])
            plt.show() ### plot to check

        totcnts = np.nansum(np.nansum(cutim))#/cutim.size
        ### Divide by itime to get flux
        A[i] = totcnts/cam.get_tint()
        # time.sleep(0.1)

    return np.mean(A)

def _acquire_fiber(fibr):
    '''
    Change to a new science fiber. If not given a valid fiber number, does nothing

    Returns: 1 for sucess, 0 if no valid fiber is given
    '''
    if fibr in [1,2,3,4]:
        track.set_goal("sf{}".format(fibr))
        # wait for tracking script to update goal
        while track.get_psf_goal()[3] != fibr:
            time.sleep(.05)

        start_t = time.time()
        # wait for a valid psf to be within .5 pixels in x and y of target
        while True:
            # calculate the dots to print based on how long we've been running
            # dots = "."*((int(time.time()) - int(start_t))%3+1) 
            # # clear line
            # sys.stdout.write(u"\u001b[2K")
            # sys.stdout.write("\rfiber {} | acquiring".format(fibr)+dots)
            # sys.stdout.flush()
            # get location of psf
            params = track.get_psf_parameters()
            # get target location
            goal = track.get_psf_goal()

            # if we're close enough to goal, break
            if params[0] and abs(params[2] - goal[1]) < .5 and abs(params[3] - goal[2]) < .5:
                return 1
    else:
        return 0

def _take_backgrounds():
    '''Take PD backgrounds

    Returns: Background PD voltage (float)
    '''
    ### Stop tracking and move off-slit for backgrounds
    TrackFlag = track.is_tracking()
    if TrackFlag:
            track.stop_tracking()
            # Wait for the Tracking loop to complete its last iteration
            time.sleep(0.5)
    ### Get initial FAM position
    FAM_CP = FAM.get_target()
    ### Move FAM and wait for new position
    FAM.set_pos([9, 9], block=True)

    ### Take backgrounds
    with redPM_cmds() as pd:
        bkgd = pd.read_pd(nreads)*1e4

    ### Return to original FAM position and restart tracking
    FAM.set_pos([FAM_CP[0], FAM_CP[1]], block=True)
    if TrackFlag:
        track.start_tracking()

    return bkgd

def _take_fam_scan(start=-3, stop=3, step=1, fibers=None, bkgd=None, Path=None):
    ''' Takes 2D FIU FAM scan

    Keyword arguments:
    start: (float, default -3) lower limit for scan in cred2 pixels
    stop:  (float, default 3) upper limit for scan in cred2 pixels
    step:  (float, default 1) stepsize of scan in cred2 pixels
    fibers: (list, default None) list of fibers to scan. If None, scans around current TT location
    bkgd: (float, default None) background value to subtract. If None, takes new background
    Path: (2-tuple string, default None) Path/additional filename to save scan results. Default is '/home/nfiudev/dev/FAM_Scans/'
    
    Outputs:
    Saves a numpy archive of the scan xpoints, ypoints, pd readings, and cred2 fluxes. Default path is 
    /home/nfiudev/dev/FAM_Scans/<timestamp>_pd_fiu_scan.npz. Use Path to override directory or add a file prefix.

    Returns:
    Tuple of xpoints, ypoints, pdfluxes, cred2fluxes. If one fiber is scanned, dimensions are 
    1, 1, 2, 2, respectively. If multiple fibers, dimensions are 1, 1, 3, 3.
    '''
    if Path is None:
        Path = ('/home/nfiudev/dev/FAM_Scans/', '')
    Track_Flag = track.is_tracking()
    ### This isn't working
    if Track_Flag:
        goal0 = int(track.get_goal()[0])
    try:
        if bkgd is None:
            bkgd = _take_backgrounds()
        if fibers is None:
            currPos=True
            # goal0 = int(track.get_goal()[0])
            fibers = [-1]
        else:
            currPos=False
        
        datacubes = []
        crd2cubes = []
        for fibr in fibers:
            ### Get on target, then turn off tracking and just use ttm
            if Track_Flag or not currPos:
                if not track.is_tracking(): track.start_tracking()
                _ = _acquire_fiber(fibr)
                track.stop_tracking()

            ### Use average position for scan center
            ipos = np.empty((2,10))
            for i in range(10):
                pos = FAM.get_pos()
                ipos[0,i] = pos[0]
                ipos[1,i] = pos[1]
                time.sleep(0.05)
            ipos = np.mean(ipos, axis=-1)

            xpixtorad = 0.01515151515
            ypixtorad = 0.01709401709
            xpts = xpixtorad*np.arange(start, stop+step/2., step)+ipos[0]
            ypts = ypixtorad*np.arange(start, stop+step/2., step)+ipos[1]
            
            FAM.set_pos([ipos[0], ipos[1]], block=True)

            datacube = np.empty((ypts.size, xpts.size))
            crd2cube = np.empty((ypts.size, xpts.size))
            ### Scan over 2D array
            cnt = 0
            nframes = str(xpts.size*ypts.size)
            with redPM_cmds() as pd:
                for i in range(ypts.size):
                    for j in range(xpts.size):
                        ### Move FIU FAM and wait to update
                        FAM.set_pos([xpts[i], ypts[j]], block=True)
                        ### I'm not convinced keeping the FIU FAM in near-constant motion doesn't cause problems
                        # time.sleep(0.05)
                        ### Read the PD
                        datacube[i,j] = 1e4*pd.read_pd(nreads)-bkgd

                        txt = 'X = %5.2f - Y = %5.2f,  Frame %3.0f / '%(float(xpts[i]), float(ypts[j]), float(cnt+1)) + nframes
                        sys.stdout.write('\rFAM Pixel Targets: ' + txt)
                        sys.stdout.flush()

                        ### Get a cred2 flux
                        crd2cube[i,j] = _get_cred2_flux()
                        cnt+=1

            cred2fluxes = np.asarray(crd2cube)
            np.savez(Path[0]+Path[1]+'_pd_fiu_scan.npz', 
                    xpts=xpts,ypts=ypts,pddata=datacube, cred2flx=crd2cube)
            
            print('')
            print('')
            FAM.set_pos([ipos[0], ipos[1]], block=True)
            datacubes.append(datacube)
            crd2cubes.append(crd2cube)

        ### Restart tracking if it was in use before
        if Track_Flag:
            track.start_tracking()
            ### Wait for the tracking script to reacquire correctly
            time.sleep(0.5)

        if len(datacubes) == 1:
            return xpts, ypts,  datacubes[0], crd2cubes[0]
        else:
            return xpts, ypts,  datacubes, crd2cubes

    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        track.set_user_offset(offset_x=0, offset_y=0)
        FAM.set_pos([ipos[0], ipos[1]], block=True)
        if Track_Flag: 
            track.set_goal(goal0)
            track.start_tracking()
        print('Tracking reset complete')
         ### Save scan even if interrupted
        np.savez(Path[0]+Path[1]+'_pd_fiu_scan.npz', 
                    xpts=xpts,ypts=ypts,pddata=datacube, cred2flx=crd2cube)
        # print(e)
        raise e

def _analyze_fam_scan(scan=None, xpts=None, ypts=None, cred2fluxes=None, plot=True, Path=None, fibers=None, savename='', fit=True, overwrite=False):
    '''Runs analysis of FIU FAM scan data

    Keyword arguments:
    scan: (ndarray, default None) PD readings. Either 2-dim (for a single scan), or 3-dim (multiple scans)
    xpoints: (ndarray, default None) xpoints of the scan. If None, uses np.arange() of appropriate size
    ypoints: (ndarray, default None) Above but ypoints
    cred2fluxes: (ndarray, default None) Same as scan, but cred2 fluxes
    plot: (bool, default True) Whether to show a plot. On Docker, make sure xhost +local:docker was run 
          first. Even if you have, this sometimes causes segfaults.
    Path: (2-tuple of str, default None) Path/name to use to load scans. Not used at present (was used for 
          loading scans from disk, which is temporarily disabled)
    fibers: (list of int, default None) List of fiber numbers corresponding to the scans
    savename: (string, default '') String to append to end of plto filename
    fit: (bool, default True) Whether to perform a Gaussian fit on the scan data
    overwrite: (bool, default False) Doesn't do anything at present, was related to whether plots/scans should be overwritten

    Returns:
    Tuple of best FAM x position, best FAM y position. If multiple scans, tuple is a list of best x positions,
    list of best y positions.
    '''

    if Path is None:
        Path = ('/home/nfiudev/dev/FAM_Scans/', '')
    if fibers is None:
        fibers = [track.get_goal()[0]]
    # ### Use a list of cubes for multiple fibers
    # crdflxs = cred2fluxes
    # if scan is not None:
    #     scan = np.asarray(scan)
    # scans = []
    # if scan is None:
    #     for i in range(len(fibers)):
    #         cname = input('Enter scan name, fiber '+str(int(fibers[i]))+' >> ').replace('\'', '').replace(' ', '')
    #         fname = Path[0]+Path[1]+cname
    #         if not os.path.exists(fname):
    #             print('Invalid scan filename!')
    #             print(fname)
    #             return -1
    #         else:
    #             data = np.load(fname)
    #             xpts = data['xpts']
    #             ypts = data['ypts']
    #             crdflxs = data['cred2flx']
    #             scans.append(data['pddata'])
    #             Path = ('', fname.split('_fiu')[0])
    # elif scan.ndim==3:
    #     scans.append(scan)
    # else:
    #     scans = scan
    # scans = np.asarray(scans)
    
    if np.ndim(scan) == 2:
        scans = np.expand_dims(scan, 0)
        cred2fluxes = np.expand_dims(cred2fluxes, 0)
    else:
        scans = np.copy(scan)
    if cred2fluxes is None:
        cred2fluxes = np.ones(scans.shape)
    if xpts is None or ypts is None:
        xpts = np.arange(scan.shape[1])
        ypts = np.arange(scan.shape[0])
        ustr = 'step'
    else: ustr = 'pos'

    bestxs, bestys = [], []
    for i, scan in enumerate(scans):
        # print(scan)
        usavename = '_fiber'+str(int(fibers[i]))+savename
        print('Fiber', str(int(fibers[i])))
        int_flux = scan[:,::-1]
        crd_flux = cred2fluxes[i][:,::-1]

        # int_flux/=crd_flux

        bestx, besty = np.unravel_index(np.argmax(int_flux), int_flux.shape)
        bestx = xpts[bestx] ### Need to flip to match with image
        besty = ypts[::-1][besty]
        
        print('Best FAM positions (raw): X = %06.2f -- Y = %06.2f -- PD = %6.3f'%(bestx,besty,np.amax(int_flux)))
        
        if not fit:
            if plot:
                extent = [xpts[0]-0.5, xpts[-1]+0.5, ypts[0]-0.5, ypts[-1]+0.5]
                plt.imshow(int_flux/np.max(int_flux), extent=extent)
                plt.colorbar()
                plt.title('FIU FAM Offset Scan')
                plt.xlabel('X pixel offset')
                plt.ylabel('Y pixel offset')
                plt.savefig(Path[0]+Path[1]+'_fam_'+usavename+'.png', bbox_inches='tight')
                plt.show()
            bestxs.append(bestx)
            bestys.append(besty)

        else:
            Flux  = int_flux
            s1, s0 = np.max(Flux - np.min(Flux)), np.min(Flux)
            Flux = (Flux - np.min(Flux))/np.max(Flux - np.min(Flux))
            # 
            Rx, Ry = np.meshgrid(xpts, ypts)
            Ry = np.flipud(Ry)
            # Computes the dimensions of the injection map.
            dim = np.shape(Flux)
            fx  = np.poly1d(np.polyfit(np.arange(dim[0]),np.mean(Rx,0),1))
            fy  = np.poly1d(np.polyfit(np.arange(dim[1]),np.mean(Ry,1),1))
            # Reshape the injection map into a vector
            Vec = Flux.ravel()
            # Prepare a list of parameters used by the curve_fit function as a 
            # starting point.     
            px,py,wx,wy = moments(Flux)
            SP  = (np.max(Flux),px, py,wx,wy,0*np.pi , np.min(Flux))
            try:
                if np.max(int_flux) < 0.1:
                    print('Peak PD reading < 0.1 V, skipping fit')
                    raise RuntimeError           
                p,u = opt.curve_fit(Gaussian2D, dim, Vec, p0 = SP,
                    bounds=[(0, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf),
                            (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)])

                fitmod = np.reshape(Gaussian2D(dim,p[0],p[1],p[2],p[3],p[4],p[5],p[6]),dim) 
                residual = Flux - fitmod
                #
                offset_x = fx(p[1])
                offset_y = fy(p[2])
                #
                # print('Gaussfit:')
                # print('Maximum: %8.6f' %(p[0]))
                print('Best FAM positions (fit): X = %06.2f -- Y = %06.2f -- PD = %6.3f' %(offset_x,offset_y, p[0]*s1+s0))
                bestxs.append(offset_x)
                bestys.append(offset_y)
            except:
                print('Fit failed! Returning raw values')
                Flux = int_flux
                fitmod = np.nan*np.zeros(Flux.shape)
                residual = np.nan*np.zeros(Flux.shape)
                bestxs.append(bestx)
                bestys.append(besty)
                offset_x = bestx
                offset_y = besty
            
            if plot:
                fig = _plot_fam_scan(Ry, Rx, Flux.T, fitmod.T, residual.T, offset_x, offset_y, Path, usavename)
                fig.show()

    if len(bestxs) == 1:
        return bestxs, bestys
    else:
        return bestxs, bestys

def _plot_fam_scan(Rx, Ry, Flux, model, residual, offset_x, offset_y, Path, savename):
    ''' Internal function for _analyze_fam_scan used to make plots. Don't mess with this,
    call _analyze_fam_scan() instead

    Seriously don't mess with this, it's taken hours to get the axis extents/orientations right!
    '''
    if np.any(np.isnan(model)):
        fig = plt.figure()
        plt.set_cmap('RdBu')
        # Add the title to this sub image
        plt.title('Injection Map')
        # Plot limits
        # Plot limits
        xstep = (np.max(Rx) - np.min(Rx))/Rx.size
        ystep = (np.max(Ry) - np.min(Ry))/Ry.size
        # print(xstep, ystep)
        limits = [np.min(Ry)-0.5*ystep,np.max(Ry)+0.5*ystep, np.min(Rx)-0.5*xstep,np.max(Rx)+0.5*xstep]
        # limits = [np.min(Rx),np.max(Rx),np.min(Ry),np.max(Ry)]
        # The injection map is display in linear scale
        plt.imshow(Flux,extent=limits, vmin=-1, vmax=1)
        # Modify the axis: one ticks every 3 pixel
        # X_ticks = np.arange(limits[0]+0.5, limits[1], 1)
        X_ticks = np.linspace(limits[0], limits[1], 5)
        plt.gca().set_xticks(X_ticks)
        # Y_ticks = np.arange(limits[2]+0.5, limits[3], 1)
        Y_ticks = np.linspace(limits[2], limits[3], 5)
        plt.gca().set_yticks(Y_ticks)
        plt.xlabel('PSF offset X direction (pixel)')
        plt.ylabel('PSF offset Y direction (pixel)')
        filename = '_raw_flux_map.pdf'
    else:
        fig = plt.figure(num = 1, figsize = (15,4))
        
        gs1 = GridSpec(1, 6)
        gs1.update(left=0.05, right=0.95, top = 0.84, bottom = 0.05, hspace=0.4, wspace=0.4)

        # Create a subplot for the image
        plt.subplot(gs1[:,:2])
        plt.set_cmap('RdBu')
        # Add the title to this sub image
        plt.title('Injection Map')
        # Plot limits
        # Plot limits
        xstep = (np.max(Rx) - np.min(Rx))/Rx.size
        ystep = (np.max(Ry) - np.min(Ry))/Ry.size
        # print(xstep, ystep)
        limits = [np.min(Ry)-0.5*ystep,np.max(Ry)+0.5*ystep, np.min(Rx)-0.5*xstep,np.max(Rx)+0.5*xstep]
        # limits = [np.min(Rx),np.max(Rx),np.min(Ry),np.max(Ry)]
        # The injection map is display in linear scale
        plt.imshow(np.abs(Flux),extent=limits, vmin=-1, vmax=1)
        # Modify the axis: one ticks every 3 pixel
        # X_ticks = np.arange(limits[0]+0.5, limits[1], 1)
        X_ticks = np.linspace(limits[0], limits[1], 5)
        plt.gca().set_xticks(X_ticks)
        # Y_ticks = np.arange(limits[2]+0.5, limits[3], 1)
        Y_ticks = np.linspace(limits[2], limits[3], 5)
        plt.gca().set_yticks(Y_ticks)
        plt.xlabel('PSF offset X direction (pixel)')
        plt.ylabel('PSF offset Y direction (pixel)')
        # Add a colorbar
        # cbar = plt.colorbar(mappable=mpl.cm.ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='Blues'))
        # Add the label
        # cbar.set_label('Normalized Flux', rotation=90)
        # Adjust the size of the tick and associated labels
        # cbar.set_clim(0,1)
        # cbar.set_ticks(np.arange(6)/5.)

        # Create a subplot for the image
        plt.subplot(gs1[:,2:4])
        plt.set_cmap('RdBu')
        # Add the title to this sub image
        plt.title('Gaussian fit')
        plt.plot([offset_x,offset_x],[offset_y,offset_y],'+r')
        # The injection map is display in linear scale
        plt.imshow(model,extent=limits, vmin=-1, vmax=1)
        # Modify the axis: one ticks every 3 pixel
        # X_ticks = np.arange(limits[0]+0.5, limits[1], 1)
        X_ticks = np.linspace(limits[0], limits[1], 5)
        plt.gca().set_xticks(X_ticks)
        # Y_ticks = np.arange(limits[2]+0.5, limits[3], 1)
        Y_ticks = np.linspace(limits[2], limits[3], 5)
        plt.gca().set_yticks(Y_ticks)
        plt.xlabel('PSF offset X direction (pixel)')
        plt.ylabel('PSF offset Y direction (pixel)')
        # Add a colorbar
        # cbar = plt.colorbar(mappable=mpl.cm.ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='Blues'))
        # Add the label
        # cbar.set_label('Normalized Flux', rotation=90)
        # Adjust the size of the tick and associated labels
        # cbar.set_clim(0,1)
        # cbar.set_ticks(np.arange(6)/5.)

        # Create a subplot for the image
        plt.subplot(gs1[:,4:])
        plt.set_cmap('RdBu')
        # Add the title to this sub image
        plt.title('Residual')
        # The injection map is display in linear scale
        plt.imshow(residual,extent=limits, vmin=-0.15, vmax=0.15)
        # Modify the axis: one ticks every 3 pixel
        # X_ticks = np.arange(limits[0]+0.5, limits[1], 1)
        X_ticks = np.linspace(limits[0], limits[1], 5)
        plt.gca().set_xticks(X_ticks)
        # Y_ticks = np.arange(limits[2]+0.5, limits[3], 1)
        Y_ticks = np.linspace(limits[2], limits[3], 5)
        plt.gca().set_yticks(Y_ticks)
        plt.xlabel('PSF offset X direction (pixel)')
        plt.ylabel('PSF offset Y direction (pixel)')
        # Add a colorbar
        # cbar = plt.colorbar(mappable=mpl.cm.ScalarMappable(norm=Normalize(vmin=-0.15, vmax=0.15), cmap='RdBu'))
        # Adjust the size of the tick and associated labels
        # cbar.set_clim(-0.1,0.1)
        # cbar.set_ticks([-0.10,-0.05,0.00,0.05,0.10])
        # cbar.set_ticklabels(['-0.10','-0.05',' 0.00',' 0.05',' 0.10'])

        filename = '_flux_map.pdf'
    
    fullname = Path[0]+Path[1]+filename
    print('Saving plot to:', fullname)
    plt.savefig(fullname, bbox_inches='tight', pad_inches=0.25, dpi=600)
    return fig
    # plt.show()

def _update_fiber_position(fnum, FAM_pos):
    ### Set to FAM position
    FAM.set_pos(FAM_pos, block=True)

    ### Get pixel values - take time average
    xvals, yvals = np.empty(40), np.empty(40)
    for i in range(40):
        psfpars = track.get_psf_parameters()
        xvals[i] = psfpars[2]
        yvals[i] = psfpars[3]
        time.sleep(0.02)
    xpix = np.round(np.mean(xvals),2)
    ypix = np.round(np.mean(yvals))

    track.set_fiber_loc((xpix, ypix), int(fnum))
    track._loadGoals()
    print('Updated fiber position:', track.goals[int(fnum)])

def fam_scan(start=-4, stop=4, step=1, fibers=None, savename='', plot=True, fit=True, overwrite=False, update_fiber=False):
    ''' Take and analyze a TT scan.

    Keyword arguments:
    start: (float, default -4) lower limit for scan in cred2 pixels
    stop:  (float, default 4) upper limit for scan in cred2 pixels
    step:  (float, default 1) stepsize of scan in cred2 pixels
    fibers: (list, default None) list of fibers to scan. If None, scans around current TT location
    savename: (string, default '') String to append to end of plto filename
    plot: (bool, default True) Whether to show a plot. On Docker, make sure xhost +local:docker was run 
          first. Even if you have, this sometimes causes segfaults.
    fit: (bool, default True) Whether to perform a Gaussian fit on the scan data
    overwrite: (bool, default False)  Doesn't do anything at present, was related to whether plots/scans should be overwritten
    update_fiber: (bool, default False) Prompts to update fiber location after fit. Currently only works for one fiber, and
                assumes the fiber to be updated is the tracking goal. 

    Returns: Tuple of best x FAM positions, best y FAM postions. If multiple scans, best positions are lists

    '''
    startTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
    tstamp = startTime.replace('-','')
    tstamp = tstamp.replace(':','')
    tstamp = tstamp.replace(' ','')
    tstamp = tstamp[2:-2]
    Path = ('/home/nfiudev/dev/FAM_Scans/',tstamp)
    ### Take backgrounds
    bkgd = _take_backgrounds()
    
    xpts, ypts, scans, cred2fluxes = _take_fam_scan(start=start, stop=stop, step=step, fibers=fibers, bkgd=bkgd, Path=Path)

    xoffset, yoffset = _analyze_fam_scan(scan=scans, xpts=xpts, ypts=ypts, cred2fluxes=cred2fluxes, 
                                plot=plot, fit=fit, Path=Path, savename=savename, fibers=fibers)

    if update_fiber:
        newx = np.round(xoffset[0],3)
        newy = np.round(yoffset[0],3)
        fnum = track.get_psf_goal()[-1]
        confirm = input('Update position of fiber '+str(int(fnum))+' to '+str(newx)+', '+str(newy)+'? [Y/n] >> ')
        if confirm in ['Y', 'Yes', 'yes']:
            _update_fiber_position(fnum, (xoffset[0], yoffset[0]))
        else:
            print('Fiber positions not updated')

    return xoffset, yoffset

def scan_mode(noll, start=-0.03, stop=0.03, step=0.01, grid_size=13, stepsize=2, refsrf='211129_601'):
    ''' Run scan over sepcified zernike mode
    
    Arguments:
    noll (int) Zernike noll to scan over

    Keyword arguments:
    start (float, default -0.03): Lower limit of aberattion axis. Units are whatever DM.py uses, IDK
    stop (float, default 0.03): Upper limit on aberration axis (inclusive)
    step (float, default 0.01): Stepsize along aberration axis
    gridsize (int, default 13): Size of TT scan do perform, in steps. Default gives a 13x13
    stepsize (float, default 2): Stepsize of the TT scan, in cred 2 pixels
    refsrf (string, default '211129_601') Name of reference DM surface to apply aberration to. If
            'current', uses whatever the DM shape is at the start of the scan (useful for iterative tuning)

    Outputs: Tuple of scan x points, y points, applied amplitudes, pd readings, and cred2 fluxes. Path is
    /home/nfiudev/dev/FAM_Scans/Noll_<noll>_PD_Values.npz. Note that this will overwrite any preivous scan of that noll!!!

    Returns: Tuple of scan x points, y points, applied amplitudes, pd readings, and cred2 fluxes.
    '''
    ### Reset to initial DM shape
    # K2AO.load_dm_shape(dmfname)
    # inimap = K2AO.DM_shape_map
    ### Loop through points
    try:
        ### Get DM reference. Either a flatmap or 'current', do do everything relative to the current
        ### DM setting
        if refsrf != 'current':
            flat = DM.setFlatSurf(refsrf)
        else:
            flat = np.copy(DM.shm_Surf.get_data())
        bkgd = _take_backgrounds()
        zpts = np.arange(start, stop+step/2., step)
        scans = np.empty((zpts.size, grid_size, grid_size))
        cred2s = np.empty((zpts.size, grid_size, grid_size))
        for i in range(zpts.size):
            print('Amplitude:', np.round(zpts[i],3))
            ### Apply zernike
            dmmap = DM.pokeZernike(zpts[i], noll, bias=flat)
            ### Read photometer
            xpts, ypts, scans[i], cred2s[i] = _take_fam_scan(bkgd=bkgd, start=-int((stepsize*(grid_size-1))/2), stop=int((stepsize*(grid_size-1))/2), step=stepsize)
            ### Set DM back where it started
            DM.setSurf(flat)

        ### Save data
        Path = '/home/nfiudev/dev/FAM_Scans/'
        fullname = Path + 'Noll_%02d_PD_Values.npz' %(noll)
        print('Saving to:', fullname)
        np.savez(fullname, xpts=xpts, ypts=ypts, amps=zpts, pdfluxes=scans, crdfluxes=cred2s)
        ### return
        return xpts, ypts, zpts, scans, cred2s

    ### Reload inital DM shape if something went wrong
    except(Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        DM.setSurf(flat)
        print('DM shape reset')
        raise e

def plot_zernike_scan_peaks(noll, fit=True):
    ''' Plots raw peak voltage of scan vs applied amplitude

    Arguments: noll (int): noll index to load. Assume default filename used by scan_mode()

    Keyword Arguments: fit (bool, default True): Whether to fit a quadratic
    '''
    Path = '/home/nfiudev/dev/FAM_Scans/'
    fullname = Path + 'Noll_%02d_PD_Values.npz' %(noll)
    scandata = np.load(fullname)
    xpts, ypts = scandata['xpts'], scandata['ypts']
    amps = scandata['amps']
    pdfluxes = scandata['pdfluxes']
    crdfluxes = scandata['crdfluxes']

    peaks = np.empty(pdfluxes.shape[0])

    for i in range(peaks.size):
        peaks[i] = np.max(pdfluxes[i])
    plt.scatter(amps, peaks, color='k')
    if fit:
        mpts = np.linspace(np.nanmin(amps), np.nanmax(amps), 100)
        gauss = models.Gaussian1D(amplitude=(np.nanmax(peaks)-np.nanmin(peaks)), mean=0, stddev=0.1)
        fitter = modeling.fitting.LevMarLSQFitter()
        fitted = fitter(gauss, amps[np.isfinite(peaks)], peaks[np.isfinite(peaks)])
        fitcurve = fitted(mpts)
        plt.plot(mpts, fitcurve, color='k', linestyle='--')
        bestpt = mpts[np.argmax(fitcurve)]
        plt.axvline(bestpt, color='r', linestyle='--')
        plt.text(bestpt+0.015, np.min(peaks)+0.5*(np.max(peaks)-np.min(peaks)), 'Fit peak = '+str(np.round(bestpt,3))+' V')
    plt.title(nolls[noll] + ' Raw Voltage')
    plt.xlabel('Mode Amplitude')
    plt.ylabel('Peak PD voltage')
    plt.savefig(Path+'noll'+str(noll)+'_rawpeaks.png', bbox_inches='tight')
    plt.show()

    if fit: return bestpt
    return amps[np.argmax(peaks)]

def plot_zernike_scans(noll):
    '''Plots 2D images of each frame in scan

    Arguments: noll (int): noll index to load. Assume default filename used by scan_mode()
    '''
    Path = '/home/nfiudev/dev/FAM_Scans/'
    fullname = Path + 'Noll_%02d_PD_Values.npz' %(noll)
    scandata = np.load(fullname)
    xpts, ypts = scandata['xpts'], scandata['ypts']
    xlow, xhigh = np.min(xpts), np.max(xpts)
    ylow, yhigh = np.min(ypts), np.max(ypts)
    bounds = [xlow, xhigh, ylow, yhigh]
    amps = scandata['amps']
    pdfluxes = scandata['pdfluxes']
    crdfluxes = scandata['crdfluxes']

    nrows = int(np.ceil(amps.size/4))

    fig, axes = plt.subplots(nrows=nrows, ncols=4)
    cnt = 0
    for row in range(nrows):
        for col in range(4):
            if cnt < amps.size:
                axes[row,col].set_title('Amp = '+str(np.round(amps[cnt],3))+' V')
                axes[row,col].imshow(np.flipud(pdfluxes[cnt]), extent=bounds)
                cnt+=1
            else:
                axes[row,col].set_axis_off()
    fig.suptitle(nolls[noll]+' Raw TT Scans')
    plt.subplots_adjust(hspace=0.01, wspace=0.35)
    plt.savefig(Path+'Noll_'+str(noll).zfill(2)+'_raw_TT_scans.png',bbox_inches='tight')
    fig.show()


def fit_scan(xpts, ypts, flux):
    ''' Fits a 2D gaussian to the given array. Expect zero flux at large distance

    Arguments:
    xpoints (ndarray): x-points of scan. Should work out to x-points on cred2 images
    ypoints (ndarray): y-points of scan
    flux (ndarray): 2D array of flux points to fit

    Returns:
    ndarray of best-fit parameters for the 2D gaussian
    '''
    ### Fit a gaussian
    Rx, Ry = np.meshgrid(xpts, ypts)
    Ry = np.flipud(Ry)
    # Computes the dimensions of the injection map.
    dim = np.shape(flux)
    fx  = np.poly1d(np.polyfit(np.arange(dim[0]),np.mean(Rx,0),1))
    fy  = np.poly1d(np.polyfit(np.arange(dim[1]),np.mean(Ry,1),1))
    # Reshape the injection map into a vector
    Vec = flux.ravel()
    # Prepare a list of parameters used by the curve_fit function as a 
    # starting point.     
    px,py,wx,wy = moments(flux)
    SP  = (np.max(flux),px, py,wx,wy,0*np.pi , np.min(flux))           
    p,u = opt.curve_fit(Gaussian2D, dim, Vec, p0 = SP)
    ### Return best-fit parameters
    return p

def plot_zernike_fit_peaks(noll, fit=True, cred2norm=False):
    '''Plots the peak of the best-fit gaussian vs applied amplitude

    Arguments: noll (int): noll index to load. Assume default filename used by scan_mode()

    Keyword Arguments: 
    fit (bool, default True): Whether to fit a quadratic to the amplitudes
    cred2norm (bool, default False): Whether to normalize by cred2 fluxes. Causes more problems lately
    '''
    Path = '/home/nfiudev/dev/FAM_Scans/'
    fullname = Path + 'Noll_%02d_PD_Values.npz' %(noll)
    scandata = np.load(fullname)
    xpts, ypts = scandata['xpts'], scandata['ypts']
    amps = scandata['amps']
    pdfluxes = scandata['pdfluxes']
    crdfluxes = scandata['crdfluxes']

    peaks = np.empty(pdfluxes.shape[0])

    for i in range(peaks.size):
        pars = fit_scan(xpts, ypts, pdfluxes[i])
        peaks[i] = pars[0]

    plt.scatter(amps, peaks, color='k')
    if fit:
        mpts = np.linspace(np.nanmin(amps), np.nanmax(amps), 100)
        gauss = models.Gaussian1D(amplitude=(np.nanmax(peaks)-np.nanmin(peaks)), mean=0, stddev=0.1)
        fitter = modeling.fitting.LevMarLSQFitter()
        fitted = fitter(gauss, amps[np.isfinite(peaks)], peaks[np.isfinite(peaks)])
        fitcurve = fitted(mpts)
        plt.plot(mpts, fitcurve, color='k', linestyle='--')
        # poly = np.poly1d(np.polyfit(amps, peaks, deg=2))
        # fitcurve = poly(mpts)
        # plt.plot(mpts, fitcurve, color='k', linestyle='--')
        bestpt = mpts[np.argmax(fitcurve)]
        plt.axvline(bestpt, color='r', linestyle='--')
        plt.text(bestpt+0.015, np.nanmin(peaks)+0.5*(np.nanmax(peaks)-np.nanmin(peaks)), 'Fit peak = '+str(np.round(bestpt,3))+' V')
    plt.title(nolls[noll]+' Best-fit Amplitude')
    plt.xlabel('Mode Amplitude')
    plt.ylabel('Best-Fit Peak PD voltage')
    plt.savefig(Path+'noll'+str(noll)+'_fitpeaks.png', bbox_inches='tight')
    plt.show()

    if fit: return bestpt
    return amps[np.argmax(peaks)]

def set_dm_shape(noll, zpt, refsrf='211129_601'):
    ''' Set DM shape to specified noll/amplitude

    Arguments:
    noll (int): noll index of zernike to apply
    zpt (float): amplitude of aberration to apply, in DM.py units

    Keyword arguments:
    refsrf (str, default '211129_601'): DM reference surface to apply
            aberration on top of. If 'current', uses the current shape
            of the DM (useful for iterative tuning)
    '''
    ### Reset to initial DM shape
    if refsrf != 'current':
        flat = DM.setFlatSurf('211129_601')
    else:
        flat = np.copy(DM.shm_Surf.get_data())
    # K2AO.load_dm_shape(dmfname)
    # inimap = K2AO.DM_shape_map
    ### Loop through points
    try:
        print('Amplitude:', np.round(zpt,3))
        ### Apply zernike
        dmmap = DM.pokeZernike(zpt, noll, bias=flat)

        return dmmap

    ### Reload inital DM shape if something went wrong
    except(Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        DM.setFlatSurf('211129_601')
        print('DM shape reset')
        raise e


def reset_dm_shape():
    '''Resets the DM shape to the 211129_601 flat map, in case you fuck up
    '''
    DM.setFlatSurf('211129_601')