#--- IMPORTS
from PD_cmds import PD_cmds
from FAM_cmds import FAM_cmds
import time
import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import warp_polar
from FIU_Fiber_cmds import FIU_Fiber_cmds
from Mode_Change_cmds import Mode_Change_cmds
from Filter_Wh_cmds import Filter_Wh_cmds
#from TCP_cmds import TCP_cmds
from PyWFS_Pickoff_cmds import PyWFS_Pickoff_cmds
from Light_Src_cmds import Light_Src_cmds
#from PIAA_cmds import PIAA_cmds
from Track_Cam_cmds import TC_cmds
from Coronagraph_cmds import Coronagraph_cmds
from SAM_cmds import SAM_cmds
from PSM_cmds import PSM_cmds
#from FM_cmds import FM_cmds
import ktl
import sys
import os
from astropy.io import fits
from pipython import GCSDevice
#sys.path.append('/home/nfiudev/dev/')
from DM_Sock import DM as dmlib
#from redPM_cmds import redPM_cmds
from serial import Serial
import atexit


plt.ion()


#--- Measurement Settings
#- PD Settings
NREAD_dflt = 10    # default number of PD reads at a given point 
VRANGE_dflt = 10   # default voltage range on PD

#--- Instantiate Devices 
TTM = FAM_cmds()
DM = dmlib()
tc = TC_cmds()
sam = SAM_cmds()
psm = PSM_cmds()
multiport = FIU_Fiber_cmds()
# ktlservice for KPIC needed for flip mirror
ktlserv = ktl.Service("kpic")


#-- Helper functions for 2micron laser comms
global las
las = None
las_usb = '/dev/serial/by-id/usb-FTDI_FT232R_USB_UART_AB0PUBJ6-if00-port0'

def laser_connect():
    global las
    las = Serial(las_usb, 115200)
    atexit.register(laser_disconnect)
    # Send message and read a couple times to synchronize
    las.write(b'power?\r')
    las.read_all()
    las.read_all()
    laser_getenable()
    las.read_all()

def laser_disconnect():
    global las
    if las is not None:
        if las.isOpen():
            las.close()
            print('2 micron laser disconnected')

def laser_query(qval):
    msg = qval + '\r'
    msg = msg.encode()
    las.write(msg)
    time.sleep(0.1)
    return las.read_all()

def laser_enable():
    print(laser_query('enable=1'))

def laser_disable():
    print(laser_query('enable=0'))

def laser_getenable():
    print(laser_query('enable?'))

def laser_getcurrent():
    print(laser_query('current?'))

def laser_setcurrent(cur):
    if cur > 450:
        raise ValueError('cannot set current above 450')
    print(laser_query('current=%i'%cur))
    
#-- Other Functions
def flip_mirror_in():
    #with FM_cmds() as flip:
    #    flip.flip(1)
    ktlserv['PDPONAM'].write('in')

def flip_mirror_out():
    #with FM_cmds() as flip:
    #    flip.flip(2)
    ktlserv['PDPONAM'].write('out')

def read_flip_mirror_state():
    #with FM_cmds() as flip:
    #    val = flip.get_state()
    return ktlserv['PDPONAM'].read()
#    return val

def read_pd(nread):
    with PD_cmds() as pd:
            pdread = pd.readN(nread)
    return pdread
    
def ttm_line_scan(ttm_dels, axis, pause, start_ttm, nread=NREAD_dflt):
    pd_reads = np.full((len(ttm_dels)),np.nan)
    for ind, ttm_del in enumerate(ttm_dels):
        if axis == 'x':
            newpos = start_ttm+np.array([ttm_del, 0])
        elif axis == 'y':
            newpos = start_ttm+np.array([0, ttm_del])
        else:
            raise ValueError("unrecognized axis '{}'".format(axis))
        TTM.set_pos(newpos)
        time.sleep(0.1)
        pd_reads[ind] = read_pd(nread).mean()
        time.sleep(pause)
    return pd_reads

def ttm_2D_scan(start, stop, nsteps, start_ttm=None, pause=0.0, nread=NREAD_dflt, VRange=VRANGE_dflt):
    '''Function to perfrom a 2D FAM scan
    
    Example: pd_reads, start_ttm, ttm_dels = tools.ttm_2D_scan(-200, 200, 20)
    '''

    # Create deltas vector
    ttm_dels = np.linspace(start, stop, nsteps)
    # Prellocate data array
    pd_reads = np.full((nsteps,nsteps),np.nan)
    # Get starting FAM position (if not provided)
    if start_ttm is None:
        start_ttm = TTM.get_pos()
    # Get X and Y positions for scan
    ttm_delxs = ttm_dels+start_ttm[0]
    ttm_delys = ttm_dels+start_ttm[1]

    # Iterate through positions, sampling PD power
    with PD_cmds(VRange=VRange) as pd:
        for xind, ttm_delx in enumerate(ttm_delxs):
            for yind, ttm_dely in enumerate(ttm_delys):
                newpos = np.array([ttm_delx, ttm_dely])
                TTM.set_pos(newpos)
                time.sleep(pause)
                pd_reads[xind,yind] = pd.readN(nread).mean()
    
    # Set FAM back to center of scan
    TTM.set_pos(start_ttm)
            
    return pd_reads, start_ttm, ttm_dels

def sam_2D_scan(start, stop, nsteps, start_sam=None, pause=0.0, nread=NREAD_dflt, VRange=VRANGE_dflt):
    # Example: pd_reads, start_sam, sam_dels = tools.sam_2D_scan(-0.2, 0.2, 20)

    sam_dels = np.linspace(start, stop, nsteps)
    
    pd_reads = np.full((nsteps,nsteps),np.nan)
    all_pos = np.full((nsteps,nsteps), 0, dtype=object)

    if start_sam is None:
        start_sam = sam.get_pos()

    sam_delxs = sam_dels+start_sam[0]
    sam_delys = sam_dels+start_sam[1]

    cur_max = -9999    
    with PD_cmds(VRange=VRange) as pd:
        for xind, sam_delx in enumerate(sam_delxs):
            #newposx = start_ttm + np.array([ttm_delx, 0])
            for yind, sam_dely in enumerate(sam_delys):
                #newpos = newposx + np.array([0, ttm_dely])
                newpos = np.array([sam_delx, sam_dely])
                print('Setting {}'.format(newpos))
                sam.set_pos(newpos, block=True)
                all_pos[xind, yind] = newpos

                time.sleep(pause)
                read_val = pd.readN(nread).mean()
                pd_reads[xind,yind] = read_val
                if read_val > cur_max:
                    cur_max = read_val
                print('---> Power = %0.5f,  (current max = %0.5f)'%(pd_reads[xind, yind], cur_max))
                
    # find SAM position of maximum power
    max_inds_1d = np.argmax(pd_reads)
    max_inds = np.unravel_index(max_inds_1d, shape=np.shape(pd_reads))
    max_pow_pos = all_pos[max_inds]

    print('Position of maximum PD power: ')
    print(max_pow_pos)
    
    # set SAM directly to max power position
    sam.set_pos(max_pow_pos, block=True)

    #return max_pow_pos
    return pd_reads, start_sam, sam_dels

def psmf_PD_scan(start, stop, nsteps, start_psmf=None, pause=0.0, nread=NREAD_dflt, VRange=VRANGE_dflt):
    # Example: pd_reads, start_psmf, psmf_dels = tools.psmf_PD_scan(-0.2, 0.2, 11)

    psmf_dels = np.linspace(start, stop, nsteps)
    
    pd_reads = np.full((nsteps),np.nan)
    
    if start_psmf is None:
        start_psmf = psm.get_foc_pos()

    psmf_dels = psmf_dels+start_psmf[0]

    with PD_cmds(VRange=VRange) as pd:
        for ind, psmf_del in enumerate(psmf_dels):
            psm.set_foc_pos(psmf_del, block=True)
            time.sleep(pause)
            pd_reads[ind] = pd.readN(nread).mean()
            
    psm.set_foc_pos(start_psmf)
            
    return pd_reads, start_psmf, psmf_dels

def multiport_PD_scan(start, stop, nsteps, start_mp=None, pause=0.0, nread=NREAD_dflt, VRange=VRANGE_dflt):
    '''Function to scan FIU_Fiber (multiport) (only 1 axis)

    Example: pd_reads, start_mp, mp_dels = tools.multiport_PD_scan(-0.2, 0.2, 11)
    '''
    mp_dels = np.linspace(start, stop, nsteps)
    
    pd_reads = np.full((nsteps),np.nan)
    
    if start_mp is None:
        start_mp = multiport.get_pos()

    mp_dels = mp_dels+start_mp[0]

    with PD_cmds(VRange=VRange) as pd:
        for ind, mp_del in enumerate(mp_dels):
            multiport.set_pos(mp_del, block=True)
            time.sleep(pause)
            pd_reads[ind] = pd.readN(nread).mean()
            
    multiport.set_pos(start_mp)
            
    return pd_reads, start_mp, mp_dels
    

def plot_2D_scan(scan, ttm_dels, start_ttm, norm=1, cropVal=7, isPlotLog=False, isCoupPlot=False):
    '''Function to plot a 2D scan

    Example: plotres = tools.plot_2D_scan(pd_reads, ttm_dels, start_ttm, cropVal=7, isCoupPlot=True)

    Returns:
    if isCoupPlot:
        return fig, ax_pix, ax_ttm , null, coords
    else:
        return fig, ax_pix, ax_ttm 
    '''
    scale_x = TTM.pix2ttmUx*4.1
    scale_y = TTM.pix2ttmUy*4.1
    
    if isPlotLog:    
        fig, (ax_pix, ax_ttm, ax_log) = plt.subplots(1,3, figsize=(17,4.91))
        plt.tight_layout()
    else:
        fig, (ax_pix, ax_ttm) = plt.subplots(1,2, figsize=(14.42,4.62))
    # .T to transpose image and get plot axes to match TTM axes
    pcm = ax_pix.imshow(scan.T/norm, origin='lower', extent=[ttm_dels[0]/scale_x, ttm_dels[-1]/scale_x, ttm_dels[0]/scale_y, ttm_dels[-1]/scale_y])
    ax_pix.set_title('KPIC 2D Coup - Center = {}'.format(start_ttm))
    ax_pix.set_ylabel('Y-Displacement [$\lambda_H/D$]')
    ax_pix.set_xlabel('X-Displacement [$\lambda_H/D$]')
    fig.colorbar(pcm, ax=ax_pix)

    pcm = ax_ttm.imshow(scan.T/norm, origin='lower', extent=[ttm_dels[0], ttm_dels[-1], ttm_dels[0], ttm_dels[-1]])
    ax_ttm.set_title('KPIC 2D Coup - Center = {}'.format(start_ttm))
    ax_ttm.set_ylabel('Y-Displacement [TTM Units]')
    ax_ttm.set_xlabel('X-Displacement [TTM Units]')
    fig.colorbar(pcm, ax=ax_ttm)

    if isPlotLog:
        pcm = ax_log.imshow(np.log10(scan.T/norm), origin='lower', extent=[ttm_dels[0], ttm_dels[-1], ttm_dels[0], ttm_dels[-1]])
        ax_log.set_title('KPIC 2D Coup (log10) - Center = {}'.format(start_ttm))
        ax_log.set_ylabel('Y-Displacement [TTM Units]')
        ax_log.set_xlabel('X-Displacement [TTM Units]')
        fig.colorbar(pcm, ax=ax_log)

    if isCoupPlot:
        null, coords = get_eta_s(scan/norm, cropVal=cropVal)
    
        ax_pix.scatter(ttm_dels[coords[0]]/scale_x, ttm_dels[coords[1]]/scale_y, s=50, c='red', marker='x')
        ax_ttm.scatter(ttm_dels[coords[0]], ttm_dels[coords[1]], s=50, c='red', marker='x')
        if isPlotLog:
            ax_log.scatter(ttm_dels[coords[0]], ttm_dels[coords[1]], s=50, c='red', marker='x')

    if isCoupPlot:
        return fig, ax_pix, ax_ttm , null, coords
    else:
        return fig, ax_pix, ax_ttm 
    
def plot_psmf_PD_scan(scan, psmf_dels, start_psmf, norm=1, cropVal=7):
    fig = plt.figure()
    pcm = plt.plot(psmf_dels, scan/norm)#extent=[ttm_dels[0]/scale_x, ttm_dels[-1]/scale_x, ttm_dels[0]/scale_y, ttm_dels[-1]/scale_y])
    plt.title('KPIC PD Focus Scan - Center = {}'.format(start_psmf))
    plt.xlabel('Focus Position [mm]')
    plt.ylabel('Power [PD Volts]')

    return fig

def plot_zernike_scan(nulls, peaks, amps, relTints=None):    
    '''Makes a simple plot of nulls and peaks vs. input Zernike amplitude

    relTints    (optional) argument with relative integration time to plot.
                    if not included, fig will only have 2 subplots
    '''
    if relTints is not None:
        fig, (ax_nulls, ax_peaks, ax_tints) = plt.subplots(1,3, figsize=(18, 4.62))
    else:
        fig, (ax_nulls, ax_peaks) = plt.subplots(1,2, figsize=(14.42, 4.62))
    
    ax_nulls.plot(amps, nulls)
    ax_nulls.set_title('KPIC Nulls vs. Zernike Scan')
    ax_nulls.set_ylabel('Null Depth [V]')
    ax_nulls.set_xlabel('Zernike Amplitude [RMS BMC units]')
    
    ax_peaks.plot(amps, peaks)
    ax_peaks.set_title('KPIC Peak Rad Avg. vs. Zernike Scan')
    ax_peaks.set_ylabel('Radially-Averaged Peak [V]')
    ax_peaks.set_xlabel('Zernike Amplitude [RMS BMC units]')

    if relTints is not None:
        ax_tints.plot(amps, relTints)
        ax_tints.set_title('Relative itime vs. Zernike Scan')
        ax_tints.set_ylabel('eta_s / eta_p^2 [1/V]')
        ax_tints.set_xlabel('Zernike Amplitude [RMS BMC units]')
        
        return fig, ax_nulls, ax_peaks, ax_tints
    else:
        return fig, ax_nulls, ax_peaks
    
def scan_zernike(noll, amp_start, amp_stop, amp_nsteps, start, stop, nsteps, start_surf=None, start_ttm=None, pause=0.0, nread=NREAD_dflt, isNorm=False):
    # Example: coup_maps, norms, amps, start_ttm, ttm_dels = tools.scan_zernike(4, -0.1, 0.1, 7, -0.2, 0.2, 20, isNorm = False)
    amps = np.linspace(amp_start, amp_stop, amp_nsteps)

    coup_maps = np.full((amp_nsteps,nsteps,nsteps), np.nan)
    norms = np.full((amp_nsteps), np.nan)
    
    if start_surf is None:
        start_surf = DM.getSurf()
    
    for ind, amp in enumerate(amps):
        print('itr: {}/{} | amp: {}'.format(ind+1, len(amps), amp))
        DM.pokeZernike(amp, noll, bias=start_surf)
        coup_map, start_ttm, ttm_dels = ttm_2D_scan(start, stop, nsteps, start_ttm=start_ttm, pause=pause, nread=nread)
        if isNorm:
            norms[ind] = meas_norm(nread)
        coup_maps[ind,:,:] = coup_map
        
    DM.setSurf(start_surf)
        
    return coup_maps, norms, amps, start_ttm, ttm_dels
    
def get_eta_s(nrmDat, cropVal=7):
    '''Function to find the null in a 2D VFN scan
        (Assumes null is within central region)
    
    Example: null, coords = tools.get_eta_s(pd_reads, cropVal=7)
    '''
    # Example for Zernike Scan check:
    # nulls = []
    # coords = []
    # for ind in range(len(amps)):
    #     null, coord = tools.get_eta_s(coup_maps[ind])
    #     nulls.append(null)
    #     coords.append(coord)
    # nulls = np.array(nulls)
    
    rowmin = np.max(int(nrmDat.shape[0]/cropVal),0)
    colmin = np.max(int(nrmDat.shape[1]/cropVal),0)
    rowmax = int((cropVal-1)*nrmDat.shape[0]/cropVal);
    colmax = int((cropVal-1)*nrmDat.shape[1]/cropVal);
    
    nrmDatCr = nrmDat[rowmin:rowmax,colmin:colmax]
    
    eta_s = np.min(nrmDatCr)
    mnInd = np.argmin(nrmDatCr)
    I1, I2 = np.unravel_index(mnInd, nrmDatCr.shape)
    
    rowMn = I1+rowmin
    colMn = I2+colmin
    
    return eta_s, (rowMn, colMn)
    
#def pol2cart(rho, phi):
#    x = rho * np.cos(phi)
#    y = rho * np.sin(phi)
#    return(x, y)
#def polarTransform(input_image, center_vec, rmax, numRadPts, numAngles):
#
#    rvec = np.linspace(0,rmax,numRadPts) # make array of radial points
#    qvec = np.linspace(0,2*np.pi-2*np.pi/numAngles,numAngles) # make array of azimuthal points
#    radialComp = np.tile(rvec, [numAngles,1]).T # 2D array of radial points 
#    angleComp = np.tile(qvec, [numRadPts, 1]) # 2D array of azimuthal points 
#    xComp, yComp = pol2cart(radialComp,angleComp) # Transform desired polar coords into cartesian coords
#    
#    # Make image coordinates 
#    rows, cols = input_image.shape
#    
#    xvals = np.arange(rows) - center_vec[0]
#    yvals = np.arange(cols) - center_vec[1]
#    
#    # compute polar transform
#    #ip = interp2d(xvals, yvals, input_image, kind='linear')
#    return ndimage.map_coordinates(input_image, [xComp.ravel(), yComp.ravel()], order=1).reshape(xComp.shape) #ip(xComp, yComp)

def polarTransform(input_image, center_vec, rmax, numRadPts, numAngles):
    rvec = np.linspace(0,rmax,numRadPts)
    qvec = np.linspace(0,2*np.pi-2*np.pi/numAngles, numAngles)
    polarim = warp_polar(image=input_image, center=center_vec, radius=rmax, output_shape=(numAngles, numRadPts))
    return polarim, rvec, qvec
    
def radAverage(img, cent=None, radpts=None, angpts=None):
    '''Function to compute radial average

    Example: radavg, rvec, polarim = tools.radAverage(pd_reads, cent=coords)
    
    '''
    # Example for Zernike Scan check:
    # radAvgs = []; rvecs = []; polarims = []
    # >>> for ind in range(len(coup_maps)):
    # ...     radavg, rvec, polarim = tools.radAverage(coup_maps[ind], cent=None)
    # ...     radAvgs.append(radavg)
    # ...     rvecs.append(rvec)
    # ...     polarims.append(polarim)
    # radAvgPeaks = np.array(radAvgs).mean(1)

    if cent is None:
        cent = np.ceil(np.array(img.shape)/2)
        
    rmax = np.min((np.array(img.shape)-np.max(cent),cent))
    
    if radpts is None:
        radpts = int(rmax)
    
    if angpts is None:
        angpts = 360
    
    polarim, rvec, qvec = polarTransform(img, cent, rmax, radpts, angpts)
    
    return np.mean(polarim, 0), rvec, polarim
    
def analyzeAroundNull(scan, norm=1, cropVal=7):
    '''Function to find the null and do a radial average automatically
    
    Example: null, coords, radprof, rvec, polarim = tools.analyzeAroundNull(pd_reads, cropVal=7)
    '''
    scan = scan/norm

    null, coords = get_eta_s(scan, cropVal=cropVal)    
    
    radprof, rvec, polarim = radAverage(scan, cent=coords)
    
    return null, coords, radprof, rvec, polarim
    
# TODO: write a function that finds center of donut automatically
    

def addVFNHeaderData(hdr=None):
    if hdr is None:
        # Make an empty fits header to start
        hdr = fits.Header()
    
    return hdr
    
def takeCRED2Image(nframes=10):
    img = tc.grab_n(nframes)[0] 
    img.header = addVFNHeaderData(img.header)
    return img
    
def saveFits(fits_nm, fits_file):
    fits_file.writeto(fits_nm)
    
def array2fits(array, header=None):
    return fits.PrimaryHDU(data=array, header=header)
    

def fullZernikeScan(noll, amp_start=-0.1, amp_stop=0.1, amp_nsteps=7, start=-0.2, stop=0.2, nsteps=20, start_surf=None, start_ttm=None, pause=0.0, nread=NREAD_dflt, isNorm=False, cropVal=5, isSave=False, savenm=None):
    # Example for running this function:
    # noll=4; result = tools.fullZernikeScan(noll, amp_start=-0.25, amp_stop=0.25, amp_nsteps=15, start=-200, stop=200, nsteps=21, start_surf=vort, start_ttm=None, pause=0.03, isNorm=False, cropVal=5, isSave=True, savenm='/nfiudata/220602/vfn_zernscan_try2_Noll%02d'%noll)

    if isSave and (savenm is None):
        raise ValueError('Provide a filename in "savenm" if isSave is True')
    
    savenm += '.fits'
    if os.path.isfile(savenm):
        raise ValueError('Provided filename exists')
        
    coup_maps, norms, amps, start_ttm, ttm_dels = scan_zernike(noll, amp_start, amp_stop, amp_nsteps, start, stop, nsteps, start_surf, start_ttm, pause, nread, isNorm)

    #-------- Find Nulls in a Zernike Scan
    nulls = []; coords = []
    for ind in range(len(amps)):
        null, coord = get_eta_s(coup_maps[ind], cropVal)
        nulls.append(null)
        coords.append(coord)
    nulls = np.array(nulls)

    #-------- Find radial average peaks in Zernike Scan
    radAvgs = []; rvecs = []; polarims = []
    for ind in range(len(coup_maps)):
        radavg, rvec, polarim = radAverage(coup_maps[ind], cent=None)
        radAvgs.append(radavg)
        rvecs.append(rvec)
        polarims.append(polarim)
    radAvgPeaks = np.array(radAvgs).max(1)
        
    #-------- Plot Zernike Scan nulls and peaks
    fig, ax_nulls, ax_peaks = plot_zernike_scan(nulls, radAvgPeaks, amps)
    fig.suptitle('Noll %d'%noll)
    
    #-------- Plot maps to check quality of scans
    for cm_ind in range(amp_nsteps):
        plotres = plot_2D_scan(coup_maps[cm_ind], ttm_dels, start_ttm, isPlotLog=True, isCoupPlot=True, cropVal=cropVal)
    
    #-------- Save Results
    fits_file = None
    if isSave:
        fits_file = array2fits(coup_maps, header=tc.grab_n(1)[0].header)
        fits_file.header.append(('','',''), end = True)
        fits_file.header.append(('','--- ZERNIKE SCAN PARAMETERS ---', ''), end = True)
        fits_file.header.append(('','',''), end = True)
        fits_file.header.append(('NollInd', noll, 'Noll Index to Scan'), end = True)
        fits_file.header.append(('AmpStart', amp_start, 'Starting Amplitude for Scan'), end = True)
        fits_file.header.append(('AmpStop', amp_stop, 'Stop Amplitude for Scan'), end = True)
        fits_file.header.append(('AmpNStep', amp_nsteps, 'Number of amplitude steps in scan'), end = True)
        fits_file.header.append(('TTMStart', start, 'Lower delta value in TTM scan'), end = True)
        fits_file.header.append(('TTMStop', stop, 'Upper delta value in TTM scan'), end = True)
        fits_file.header.append(('TTMNStep', nsteps, 'Number of steps in TTM scan'), end = True)
        fits_file.header.append(('StrtTTMX', start_ttm[0], 'Center X-value in TTM scan'), end = True)
        fits_file.header.append(('StrtTTMY', start_ttm[1], 'Center Y-value in TTM scan'), end = True)
        fits_file.header.append(('Pause', pause, 'Pause time between sample points'), end = True)
        fits_file.header.append(('isNorm', isNorm, 'Flag to take normalizations'), end = True)
        fits_file.header = addVFNHeaderData(fits_file.header)
        saveFits(savenm, fits_file)
        if isNorm:
            norm_file = array2fits(norms, header=fits_file.header)
            savenm = os.path.splitext(savenm)[0]
            saveFits(savenm+'_norms.fits', norm_file)
    
    return (coup_maps, norms, amps, start_ttm, ttm_dels, nulls, coords, radAvgs, rvecs, polarims, radAvgPeaks, fig, ax_nulls, ax_peaks, fits_file)

def fullHighQualityScan(start=-200, stop=200, nsteps=71, start_ttm=None, pause=0.0, nread=NREAD_dflt, isNorm=False, cropVal=7, isSave=False, savenm=None):
    # Example: result = tools.fullHighQualityScan(start=-200, stop=200, nsteps=71, start_ttm=None, pause=0.03, isSave=True, isNorm=False, cropVal=5, savenm='/nfiudata/KPIC_VFN/20220127_TST1/KeckPup_ch1_DedicatedTunedAgainDMMap_PIAA')
    if isSave and (savenm is None):
        raise ValueError('Provide a filename in "savenm" if isSave is True')
        
    pd_reads_raw, start_ttm, ttm_dels = ttm_2D_scan(start, stop, nsteps, start_ttm, pause, nread)
    if isNorm:
        norm_raw = meas_norm()
    else:
        norm_raw = 1.

    TTM.set_pos('background')
    time.sleep(0.5)
    bkgd = read_pd(nread).mean()
    TTM.set_pos(start_ttm)

    pd_reads = pd_reads_raw.copy() - bkgd
    if isNorm:
        norm = norm_raw.copy() - bkgd
    else:
        norm = 1.

    plotres = plot_2D_scan(pd_reads, ttm_dels, start_ttm, norm=norm,cropVal=cropVal, isPlotLog=True, isCoupPlot=True)
    radavg, rvec, polarim = radAverage(pd_reads/norm)
    stepsz = ttm_dels[1]-ttm_dels[0]
    plt.figure()
    plt.plot(rvec*stepsz/TTM.pix2ttmUx/4.1, radavg*100)
    plt.xlabel('Displacement [$\lambda_H/D$]')
    plt.ylabel('Coupling [%]')
    plt.title('KPIC VFN Donut Radial Avg. Profile\nNull=%0.2e, Peak=%0.2f%%'%(plotres[-2],radavg.max()*100))   
    
    #-------- Save Results
    fits_file = None
    if isSave:
        fits_file = array2fits(pd_reads_raw, header=tc.grab_n(1)[0].header)
        fits_file.header.append(('','',''), end = True)
        fits_file.header.append(('','--- VFN SCAN PARAMETERS ---', ''), end = True)
        fits_file.header.append(('','',''), end = True)
        fits_file.header.append(('TTMStart', start, 'Lower delta value in TTM scan'), end = True)
        fits_file.header.append(('TTMStop', stop, 'Upper delta value in TTM scan'), end = True)
        fits_file.header.append(('TTMNStep', nsteps, 'Number of steps in TTM scan'), end = True)
        fits_file.header.append(('StrtTTMX', start_ttm[0], 'Center X-value in TTM scan'), end = True)
        fits_file.header.append(('StrtTTMY', start_ttm[1], 'Center Y-value in TTM scan'), end = True)
        fits_file.header.append(('StepSize', stepsz, 'Step Size for TTM Scan'), end = True)
        fits_file.header.append(('Pause', pause, 'Pause time between sample points'), end = True)
        fits_file.header.append(('MeasDev', 'PD', 'Device used for measurements (PD or redPM)'), end = True)
        fits_file.header.append(('','',''), end = True)
        fits_file.header.append(('','--- Results ---',''), end=True)
        fits_file.header.append(('','',''), end = True)
        fits_file.header.append(('isNorm', isNorm, 'Was normalization taken'), end=True)
        fits_file.header.append(('Norm', norm_raw, 'Raw Normalization value'), end = True)
        fits_file.header.append(('bkgd', bkgd, 'Background measurement'), end = True)
        fits_file.header.append(('','NOTE: scan values are raw (no bkgd)',''), end = True)
        fits_file.header = addVFNHeaderData(fits_file.header)
        saveFits(savenm+'.fits', fits_file)
        
    return (pd_reads_raw, start_ttm, ttm_dels, norm_raw, bkgd, pd_reads, norm, plotres, radavg, rvec, polarim, stepsz, fits_file)
   
'''

####################### ZERNIKE SCANS ##################
#-------- Scan a single Zernike
coup_maps, norms, amps, start_ttm, ttm_dels = tools.scan_zernike(4, -0.1, 0.1, 7, -0.2, 0.2, 20, isNorm = False)

#-------- Find Nulls in a Zernike Scan
nulls = []; coords = []
for ind in range(len(amps)):
    null, coord = tools.get_eta_s(coup_maps[ind], cropVal = 5)
    nulls.append(null)
    coords.append(coord)

    
nulls = np.array(nulls)



#-------- Find radial average peaks in Zernike Scan
radAvgs = []; rvecs = []; polarims = []
for ind in range(len(coup_maps)):
    radavg, rvec, polarim = tools.radAverage(coup_maps[ind], cent=None)
    radAvgs.append(radavg)
    rvecs.append(rvec)
    polarims.append(polarim)

    
radAvgPeaks = np.array(radAvgs).max(1)
    
#-------- Plot Zernike Scan nulls and peaks
tools.plot_zernike_scan(nulls, radAvgPeaks, amps)









####################### Radially-Averaged Line Profile ###
pd_reads, start_ttm, ttm_dels = tools.ttm_2D_scan(-0.2, 0.2, 71)
norm = tools.meas_norm()

tools.TTM.set_pos('corner')
time.sleep(0.5)
bkgd = tools.read_pd(tools.NREAD_dflt).mean()

pd_reads_raw = pd_reads.copy()
norm_raw = norm.copy()

pd_reads -= bkgd
norm -= bkgd
plotres = tools.plot_2D_scan(pd_reads, ttm_dels, start_ttm, norm=norm)
radavg, rvec, polarim = tools.radAverage(pd_reads/norm)
stepsz = ttm_dels[1]-ttm_dels[0]
plt.figure()
plt.plot(rvec*stepsz/tools.TTM.pix2ttmUx/4.1, radavg*100)
plt.xlabel('Displacement [$\lambda_H/D$]')
plt.ylabel('Coupling [%]')
plt.title('KPIC VFN Donut Radial Avg. Profile\nNull=%0.2e, Peak=%0.2f%%'%(plotres[-2],radavg.max()*100))


####################### Automate multiple zernikes scanned
result = []
nolls = [4, 5, 6, 7, 8, 9, 10, 11]
for noll in nolls:
    result.append(tools.fullZernikeScan(noll=noll,amp_nsteps=11, start=-0.15, stop=0.15, start_ttm=tools.TTM_VortPos, isSave=True,savenm='/nfiudata/KPIC_VFN/20220127_TST1/zernScan_keckPup_ch1_Noll%d'%noll))


####################### Replot results from fullHighRes...
plotres = tools.plot_2D_scan(result[5], result[2],result[1],norm=result[6],cropVal=3)
radavg, rvec, polarim = tools.radAverage(result[5]/result[6], cent=plotres[-1])
tools.plt.figure()
tools.plt.plot(rvec*result[11]/tools.TTM.pix2ttmUx/4.1, radavg*100)
tools.plt.xlabel('Displacement[$\lambda_H/D$]')
tools.plt.ylabel('Coupling [%]')
tools.plt.title('KPIC VFN Donut Radial Avg. Profile\nNull=%0.2e, Peak=%0.2f%%'%(plotres[-2],radavg.max()*100))

'''
