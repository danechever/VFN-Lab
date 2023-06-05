'''
This library provides a function for scannling Zernikes on the DM using the PD as the sensor
'''
import numpy as np
import matplotlib.pyplot as plt 
import os, sys
from Star_Tracker_cmds import Tracking_cmds
from PD_cmds import PD_cmds
from SAM_cmds import SAM_cmds
from FAM_cmds import FAM_cmds
from DM_Sock import DM as dmlib
import time
from datetime import datetime as dt
import VFN_Toolkit_summit as tools

# Instantiate control objects
track = Tracking_cmds()
DM = dmlib()
SAM = SAM_cmds()
FAM = FAM_cmds()

# Get current date (as string for directory name)
date = dt.utcnow().strftime('%y%m%d')

############## Define Helper functions
def printv(verbose, msg):
    '''Helper function to only print if verbose=True'''
    if not verbose:
        return
    else:
        print(msg)

def DM_setup(refsrf, verbose):
    '''Helper function to set the correct DM surf'''
    if refsrf != 'current':
        printv(verbose, 'Setting custom flatmap')
        flat = DM.setFlatSurf(refsrf)
    else:
        printv(verbose, 'Using current DM surface')
        flat = DM.getSurf()
    
    return flat

def SAM_setup(fiber, verbose):
    '''Helper function to point the SAM at the right fiber on the PD'''
    # Get current SAM position
    SAM_posx, SAM_posy = SAM.get_pos()
    # Get goal position for the given fiber
    SAM_goalx, SAM_goaly = SAM.presets['pd_sf%d'%fiber]
    # Check that we are within a given tolerance (0.03)
    if ( abs(SAM_goalx - SAM_posx) > 0.03 ) or ( abs(SAM_goaly - SAM_posy) > 0.03):
        printv(verbose, 'SAM is currently at: {}, {}'.format(SAM_posx, SAM_posy))
        printv(verbose, 'Moving SAM to: {}, {}'.format(SAM_goalx, SAM_goaly))
        SAM.set_pos('pd_sf%d'%fiber, block=True)
    else:
        printv(verbose, 'SAM already in position at: {}, {}'.format(SAM_goalx, SAM_goaly))

def take_pd_background(nreads, VRange):
    '''Helper function to take a background reading on the PD
    
    NOTE: this leaves the tracking loop off but the PSF close to the starting position
    '''
    # Stop tracking and move off detector
    track.stop_tracking()
    # Wait for the Tracking loop to complete its last iteration
    time.sleep(0.5)
    # Get initial FAM position (so we can return there afterwards)
    FAM_pos0 = FAM.get_pos()
    # Move FAM
    FAM.set_pos('background', block=True)
    # wait an instant for PD voltage to settle
    time.sleep(0.2)
    # take background
    with PD_cmds() as pd:
        bkgd = pd.read_pd(nreads, v_range=VRange)   # This version returns an average already
    # Return to original FAM position to give good initial condition for tracking
    FAM.set_pos(FAM_pos0, block=True)

    return bkgd

def acquire_fiber_new(fiber, tolerance = 0.3, verbose=False):
    '''New, stripped-down function to acquire the fiber 

    NOTE: leaves tracking on at end
    
    fiber       (int) science fiber to log onto
    tolerance   (float) error tolerance to define "locked on", in pixels
    verbose     (bool) flag whether to print updates on status or keep everything quiet
    '''
    if fiber not in [1, 2, 3, 4]:
        raise ValueError("'fiber' argument must be an int identified to a science fiber")

    # Start tracking
    track.start_tracking()
    # time.sleep(1)
    # Update the fiber goal (in case it wasn't correct already)
    track.set_goal("scf{}".format(fiber))
    # wait for tracking script to update goal
    while track.get_psf_goal()[3] != fiber:
        time.sleep(0.1)
    # Wait for PSF to be within tolerance of goal
    start_t = time.time()
    while True:
        # format a live print to let users know system is still trying to lock on
        if verbose:
            dots = "."*((int(time.time()) - int(start_t))%3+1)
            sys.stdout.write(u"\u001b[2K")
            sys.stdout.write("\rfiber {} | acquiring".format(fiber)+dots)
            sys.stdout.flush()
        # Get current PSF location
        valid, _, _, psfx, psfy, _, _ = track.get_psf_parameters()
        # Get target location
        _, goalx, goaly, _ = track.get_psf_goal()
        # If we're close enough, break since we're locked
        if valid and ( abs(psfx-goalx) < tolerance ) and ( abs(psfy - goaly) < tolerance):
            break

############### Main Function
def vfn_diode_scan(noll, fiber, start=-0.3, stop=0.3, step=0.01, refsrf = 'current', nreads=10, VRange = 10, date = date, verbose=True, isSave=True, isTrackBetween=True):
    '''
    noll        Zernike to scan/tune
    fiber       (int) science fiber to use
    start       lower amplitude for noll scan
    stop        upper amplitude for noll scan
    step        amplitude step size for noll scan
    refsrf      name of flatmap to use as starting point. 'current' to use current DM surface
                  (if not 'current', must be a key in DM.flatoptdict)    
    nreads      number of PD samples to take at each point
    VRange      voltage range setting for PD (.1, 1, or 10)
    date        (str) name of subdirectory into which data will be saved - nominally a string date
    verbose     (bool) flag whether to print updates on status or keep everything quiet
    isSave      (bool) flag whether to save the output maps

    DEBUGGING ELEMENTS (for use during testing):
    isTrackBetween  (bool) flag whether tracking should be used to recenter PSF between samples

    returns: zpts, PD_volts, PD_volt0, bkgd, minIndex, minAmp
    '''

    #-- Set the DM surface
    flat = DM_setup(refsrf=refsrf, verbose=verbose)

    #-- Make sure SAM is aligned to the goal fiber
    SAM_setup(fiber=fiber, verbose=verbose)

    #-- Take PD background
    printv(verbose, 'Taking PD background')
    bkgd = take_pd_background(nreads=nreads, VRange=VRange)
    time.sleep(0.5)
    
    #-- Lock on fiber
    printv(verbose, 'Locking on Fiber')
    acquire_fiber_new(fiber, verbose=verbose)
    printv(verbose, 'Locked onto fiber %d'%fiber)

    
    #-- Save starting PD value (to get a sense of initial condition)
    track.stop_tracking()
    # wait an instant for PD voltage to finish settling
    time.sleep(0.2)
    with PD_cmds() as pd:
        PD_volt0 = pd.read_pd(nreads, v_range=VRange)   # This version returns an average already


    #-- Preallocate data arrays
    zpts  = np.round(np.arange(start, stop+step/2., step), 4)   # round to net get ridiculously small values
    PD_volts = np.zeros(len(zpts)) *np.nan     # set nonsensical value (so we know if an error occurs)


    #-- Loop through amplitudes to scan this zernike
    for ind, zpt in enumerate(zpts):
        printv(verbose, '\tAmplitude: %0.3f'%zpt)
        # Apply zernike
        DM.pokeZernike(zpt, noll, bias=flat)
        # Read PD power
        with PD_cmds() as pd:
            PD_volts[ind] = pd.read_pd(nreads, v_range=VRange)
        if isTrackBetween:
            # Set DM back to where it started and re-acquire fiber so we have a good starting condition
            DM.setSurf(flat)
            acquire_fiber_new(fiber, verbose=verbose)
            # stop tracking again for next iteration
            track.stop_tracking()

    #-- Reset the flat surface
    DM.setSurf(flat)

    #-- Analyze results
    # Subtract background
    PD_volt_clean   = PD_volts - bkgd
    PD_volt0_clean  = PD_volt0 - bkgd
    # Find where min occurs
    minIndex = np.argmin(PD_volt_clean)
    minAmp = zpts[minIndex]
    print('Minimum Power at: %0.2f'%minAmp)
    # Perform quadratic fit to nulls, focusing on +/-0.1 BMC Unit region around null only 
    indices = np.where(np.logical_and( (-0.1 <= zpts) , (zpts <= 0.1) ))
    coeffs = np.polyfit(zpts[indices], PD_volt_clean[indices], 2)

    # Plot
    plt.figure()
    plt.scatter(zpts,PD_volt_clean)
    plt.plot(zpts,[coeffs[0]*zpt**2+coeffs[1]*zpt+coeffs[2] for zpt in zpts], '--')
    plt.axvline(x=minAmp)
    plt.xlabel('Amplitude [RMS BMC Units]')
    plt.ylabel('PD Power [V]')
    titlestr = 'Zernike Scan (Noll=%d)\n'%noll
    titlestr += 'TrackBtwn = {}'.format(isTrackBetween)
    plt.title(titlestr)

    #-- Save results
    if isSave:
        # Make sure directory exists
        Path = '/nfiudata/VFN_Cals/'
        date_dir = os.path.join(Path, date)
        if not os.path.exists(date_dir):
            os.makedirs(date_dir)
        # Save results
        savenm = date_dir+'/Noll%2d_Fiber%d_time%d.png'%(noll, fiber, time.time())
        plt.savefig(savenm)
        printv(verbose, 'Saved: '+savenm)

    return zpts, PD_volts, PD_volt0, bkgd, minIndex, minAmp


########### Function that does 2D (TT) Scan 
#           ---> I'm thinking this would be a good way to start the scans

def vfn_3D_PD_scan(noll, fiber, start=-0.2, stop=0.2, step=0.05, FAM_stepRad = 200, FAM_Nstep = 20, refsrf = 'current', nreads=100, VRange = 10, delay=0.01, cropVal=7, date = date, verbose=True, isSave=True):
    '''
    noll        Zernike to scan/tune
    fiber       (int) science fiber to use
    start       lower amplitude for noll scan
    stop        upper amplitude for noll scan
    step        amplitude step size for noll scan
    FAM_stepRad FAM scan will be +/- this value 
    FAM_Nstep   Number of steps to make in FAM scan
    refsrf      name of flatmap to use as starting point. 'current' to use current DM surface
                  (if not 'current', must be a key in DM.flatoptdict)    
    nreads      number of PD samples to take at each point
    VRange      voltage range setting for PD (.1, 1, or 10)
    date        (str) name of subdirectory into which data will be saved - nominally a string date
    verbose     (bool) flag whether to print updates on status or keep everything quiet
    isSave      (bool) flag whether to save the output maps

    returns: zpts, PD_volts, PD_volt0 , bkgd, minIndex, minAmp, nulls, coords, radAvgs, start_ttm, ttm_dels
    '''

    #-- Set the DM surface
    flat = DM_setup(refsrf=refsrf, verbose=verbose)

    #-- Make sure SAM is aligned to the goal fiber
    SAM_setup(fiber=fiber, verbose=verbose)

    #-- Take PD background
    printv(verbose, 'Taking PD background')
    bkgd = take_pd_background(nreads=nreads, VRange=VRange)
    
    #-- Lock on fiber
    printv(verbose, 'Locking on Fiber')
    acquire_fiber_new(fiber, verbose=verbose)
    printv(verbose, 'Locked onto fiber %d'%fiber)

    
    #-- Save starting PD value (to get a sense of initial condition)
    track.stop_tracking()
    # wait an instant for PD voltage to finish settling
    time.sleep(0.2)
    with PD_cmds() as pd:
        PD_volt0 = pd.read_pd(nreads, v_range=VRange)   # This version returns an average already


    #-- Preallocate data arrays
    zpts  = np.round(np.arange(start, stop+step/2., step), 4)   # round to net get ridiculously small values
    PD_volts = np.zeros((len(zpts), FAM_Nstep, FAM_Nstep)) *np.nan     # set nonsensical value (so we know if an error occurs)
    # preset start_ttm to None. First iteration will then use current FAM pos and all others will use that same pos
    start_ttm = None

    #-- Loop through amplitudes to scan this zernike
    for ind, zpt in enumerate(zpts):
        printv(verbose, '\tAmplitude: %0.3f'%zpt)
        # Apply zernike
        DM.pokeZernike(zpt, noll, bias=flat)
        # Do 2D FAM scan
            # NOTE: all scans will share the same start_ttm and ttm_dels
        pd_reads, start_ttm, ttm_dels = tools.ttm_2D_scan(FAM_stepRad, -FAM_stepRad, FAM_Nstep, 
                start_ttm=start_ttm, pause=delay, nread=nreads, VRange=VRange)
        # Save results into main array
        PD_volts[ind,:,:] = pd_reads.copy()

        # NOTE: don't re-acquire fiber so that all scans are centered equally

    #-- Reset the flat surface
    DM.setSurf(flat)


    #-- Analyze results
    printv(verbose, 'Analyzing Results')
    # Subtract background
    PD_volt_clean   = PD_volts - bkgd
    PD_volt0_clean  = PD_volt0 - bkgd
    # Find nulls and do radial averaging
    nulls = []; coords = []; radAvgs = []; rvecs = []; polarims = []
    for ind in range(len(zpts)):
        null, coord, radAvg, rvec, polarim = tools.analyzeAroundNull(PD_volt_clean[ind], cropVal=cropVal)
        nulls.append(null)
        coords.append(coord)
        radAvgs.append(radAvg)
        rvecs.append(rvec)
        polarims.append(polarim)
    nulls = np.array(nulls)
    radAvgPeaks = np.array([rav.max() for rav in radAvgs])
    
    relTints = nulls / radAvgPeaks**2

    # Find where min occurs
    minIndex = nulls.argmin()
    minAmp = zpts[minIndex]   
    print('Best Null at: %0.2f'%minAmp)

    # Perform quadratic fit to nulls, focusing on +/-0.1 BMC Unit region around null only 
#    indices = [(minIndex - 0.1/step),(minIndex + 0.1/step)]
#    for i in range(len(indices)):
#        if indices[i] < 0 or indices[i] > len(zpts):
#            indices.pop(indices[i])
#    coeffs = np.polyfit(zpts[indices[0]:indices[-1]],nulls[indices[0]:indices[-1]],2)
    indices = np.where(np.logical_and( (-0.1 <= zpts) , (zpts <= 0.1) ))
    coeffs = np.polyfit(zpts[indices], nulls[indices], 2)
#    coeffs = np.polyfit(zpts,nulls,2)
  
    #-- Plots
    printv(verbose, 'Plotting Nulls and Peaks')
    # Plot the nulls and peaks (overlayed in 1 plot)
    fig_keypoints, ax_nulls, ax_peaks, ax_tints = tools.plot_zernike_scan(nulls, radAvgPeaks, zpts, relTints)
    # Overlay the curve fit on the null plot
    ax_nulls.plot(zpts,[coeffs[0]*zpt**2+coeffs[1]*zpt+coeffs[2] for zpt in zpts], '--')
    # Add vertical marker for where best null occurs on both plots
    ax_nulls.axvline(x=minAmp)
    ax_peaks.axvline(x=minAmp)
    ax_tints.axvline(x=minAmp)
    fig_keypoints.suptitle('Zernike Scan (Noll=%d)'%noll)

    # Plot 2D maps (only the two edges and middle amplitude samples)
    printv(verbose, 'Plotting 2D maps')
    fig_2Dmaps, axs = plt.subplots(1,3, figsize=(17, 4.91))
    data_inds = [0, len(zpts)//2, -1]    # indices for data maps to plot in each subplot (first, mid, last)
    for ind, ax in enumerate(axs):
        data_ind = data_inds[ind]
        # .T to transpose image and get plot axes to match TTM axes
        pcm = ax.imshow(PD_volt_clean[data_ind].T, origin='lower', extent=[ttm_dels[0], ttm_dels[-1], ttm_dels[0], ttm_dels[-1]])
        # Add null-point marker
        ax.scatter(ttm_dels[coords[data_ind][0]], ttm_dels[coords[data_ind][1]], s=50, c='red', marker='x')
        # Format plot
        ax.set_title('Amp %0.4f'%zpts[data_ind])
        ax.set_xlabel('$\Delta$X [FAM units]')
        ax.set_ylabel('$\Delta$Y [FAM units]')
        fig_2Dmaps.colorbar(pcm, ax=ax)
    # Format figure
    plt.suptitle('2D Maps (Noll = %d)'%noll)
    plt.tight_layout()

    # Plot the point with the best null
    printv(verbose, 'Plotting 2D map with best null')
    fig_2DBest = plt.figure()
    plt.imshow(PD_volt_clean[minIndex].T, origin='lower', extent=[ttm_dels[0], ttm_dels[-1], ttm_dels[0], ttm_dels[-1]])
    # Add null-point marker
    plt.scatter(ttm_dels[coords[minIndex][0]], ttm_dels[coords[minIndex][1]], s=50, c='red', marker='x')
    # Format plot
    plt.title('Best Null (Noll=%d)\nAmp %0.4f'%(noll, zpts[minIndex]))
    plt.xlabel('$\Delta$X [FAM units]')
    plt.ylabel('$\Delta$Y [FAM units]')
    plt.colorbar()

    # Plot the radial profiles
    printv(verbose, 'Plotting radial profiles')
    fig_rads = plt.figure()
    fam_stepsz = ttm_dels[1] - ttm_dels[0]
    for ind in range(len(zpts)):
        plt.plot(rvecs[ind]*fam_stepsz, radAvgs[ind], ':', label='Amp %0.4f'%zpts[ind])
    plt.legend()
    plt.xlabel('Separation [FAM Units]')
    plt.ylabel('PD Power [V]')
    plt.title('Radial Profiles (Noll=%d)'%noll)

    #-- Save results
    if isSave:
        # Make sure directory exists
        Path = '/nfiudata/VFN_Cals/'
        date_dir = os.path.join(Path, date)
        if not os.path.exists(date_dir):
            os.makedirs(date_dir)

        svtime = time.time()
        # Save plots
        savenm = date_dir+'/time%d_Noll%2d_Fiber%d_keypoints.png'%(svtime, noll, fiber)
        plt.figure(fig_keypoints)
        plt.savefig(savenm)
        printv(verbose, 'Saved: '+savenm)

        savenm = date_dir+'/time%d_Noll%2d_Fiber%d_2Dmaps.png'%(svtime, noll, fiber)
        plt.figure(fig_2Dmaps)
        plt.savefig(savenm)
        printv(verbose, 'Saved: '+savenm)

        savenm = date_dir+'/time%d_Noll%2d_Fiber%d_2DBest.png'%(svtime, noll, fiber)
        plt.figure(fig_2DBest)
        plt.savefig(savenm)
        printv(verbose, 'Saved: '+savenm)

        savenm = date_dir+'/time%d_Noll%2d_Fiber%d_radProf.png'%(svtime, noll, fiber)
        plt.figure(fig_rads)
        plt.savefig(savenm)
        printv(verbose, 'Saved: '+savenm)

    return zpts, PD_volts, PD_volt0 , bkgd, minIndex, minAmp, nulls, coords, radAvgs, start_ttm, ttm_dels

########### Function that does a single 2D scan, with analysis
def vfn_2D_validate(FAM_stepRad=200, FAM_Nstep=20, nreads=100, VRange=10, cropVal=5, delay=0.01, verbose=True):
    '''Function to do a simple 2D scan at current location and with current DM map
        and do analysis for the resulting scan

    NOTE: does NOT lock on fiber, scans around starting point
    
    FAM_stepRad FAM scan will be +/- this value 
    FAM_Nstep   Number of steps to make in FAM scan
    nreads      number of PD samples to take at each point
    VRange      voltage range setting for PD (.1, 1, or 10)
    cropVal     subwindow for finding null in frame
    delay       time to wait at each datapoint before making a read (PD settle)
    verbose     (bool) flag whether to print updates on status or keep everything quiet
    '''

    #-- Take PD background
    printv(verbose, 'Taking PD background')
    bkgd = take_pd_background(nreads=nreads, VRange=VRange)
    time.sleep(0.5)
    
    #-- Do 2D FAM scan
    pd_reads, start_ttm, ttm_dels = tools.ttm_2D_scan(FAM_stepRad, -FAM_stepRad, FAM_Nstep, 
                pause=delay, nread=nreads, VRange=VRange)

    #-- Set FAM back to center position
    FAM.set_pos(start_ttm) 

    #-- Analyze the results
    # Background subtract
    pd_reads_clean = pd_reads - bkgd
    # Find null and radial average
    null, coords, radAvg, rvec, polarim = tools.analyzeAroundNull(pd_reads_clean, cropVal=cropVal)

    #-- Plot the scan
    plt.figure()
    plt.imshow(pd_reads_clean.T, origin='lower', extent=[ttm_dels[0], ttm_dels[-1], ttm_dels[0], ttm_dels[-1]])
    plt.scatter(ttm_dels[coords[0]], ttm_dels[coords[1]], s=50, c='red', marker='x')
    plt.title('2D VFN Scan - Center = {}\nNull = {} ({})'.format(start_ttm, null, coords))
    plt.xlabel('$\Delta$X [FAM units]')
    plt.ylabel('$\Delta$Y [FAM units]')
    plt.colorbar()

    #-- Plot the radial profile
    plt.figure()
    stepsz = ttm_dels[0] - ttm_dels[1]
    plt.plot(rvec*stepsz, radAvg, 'o-')
    plt.xlabel('Separation [FAM Units]')
    plt.ylabel('PD Power [V]')
    plt.title('Radial Profile - Peak = {}'%radAvg.max())

    return pd_reads, start_ttm, ttm_dels, bkgd, null, coords, radAvg, rvec
