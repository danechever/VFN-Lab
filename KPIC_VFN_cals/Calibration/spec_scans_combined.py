print('\n-- RUNNING VERSION FROM SVN --\n\n')

from tracemalloc import start
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
# from kpicdrp import background, trace, extraction
import kpic_background as background 
import kpic_trace as trace 
import kpic_extraction as extraction
import os
from scipy import signal
import scipy.optimize as opt
import scipy
import multiprocessing as mp
from Star_Tracker_cmds import Tracking_cmds
from Track_Cam_cmds import TC_cmds
from Track_Cam_process import TC_process
from FAM_cmds import FAM_cmds
from SAM_cmds import SAM_cmds
from PSM_cmds import PSM_cmds
from Coronagraph_cmds import Coronagraph_cmds
from PyWFS_Pickoff_cmds import PyWFS_Pickoff_cmds
from FIU_Fiber_cmds import FIU_Fiber_cmds
from PyWFS_Filter_cmds import PyWFS_Filter_cmds
from Filter_Wh_cmds import Filter_Wh_cmds
from PIAA_cmds import PIAA_cmds
from Telescope_Simulator_cmds import Telescope_Simulator_cmds
from Mode_Change_cmds import Mode_Change_cmds
import Acquisition
import ktl
import Organization
from PSF_finder import moments, Gaussian2D, PSF_Finder
import sys
import time
import warnings
import datetime
import pdb
import subprocess as sp
import glob
from scipy.signal import medfilt as medfilt
# sys.path.append('/home/nfiudev/dev/')
from PD_cmds import PD_cmds as redPM_cmds
# sys.path.append('/home/nfiudev/dev/dechever/DMCharac/Code/')
from DM_Sock import DM
import VFN_Toolkit_summit as tools
from extract_spec_flux import _total_flux_frame, _calc_trace

# useful functions when observing
from Acquisition import fiber_bounce
from Scans import PS_Scan

# the light src retractor library has the same name as Light_Src_cmds used for k2ao source...
# import in this order for now
from Light_Src_cmds import Light_Src_cmds
LightSrcRet = Light_Src_cmds()
from Nirspec_cmds import Spec_cmds, Scam_cmds, Light_Src_cmds

sys.path.append('/nfiudata/Python_Libraries/')
# Library for moving SFP position
# sys.path.append('/nfiudata/Python_Libraries/SFP_cmds.py')
# sys.path.append('/nfiudata/Python_Libraries/Nirc2_cmds.py')
# TODO: Filter out AO server output when importing this library
from SFP_cmds import SFP_cmds

# Library for controlling nirc2
from Nirc2_cmds import Nirc2_cmds

'''
Required dependencies;
    kpic_background.py
    kpic_trace.py
    kpic_extraction.py
    hip81497_200602_sf2_wvsoln.txt
    2massK.txt
'''

### Get tons of these from extraction/rectification due to bad pixels
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Filter out astropy warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)

print('Current save directory:', Organization.get_path()[0])

# instantiate the stages etc.
cam   = TC_cmds()
track = Tracking_cmds()
proc = TC_process()
FAM  = FAM_cmds()
spec  = Spec_cmds()
scam  = Scam_cmds()
PSM   = PSM_cmds()
SAM   = SAM_cmds()
PIAA = PIAA_cmds()
Coronagraph = Coronagraph_cmds()
PyWFSPickoff = PyWFS_Pickoff_cmds()
PyWFSFilt = PyWFS_Filter_cmds()
FilterWh = Filter_Wh_cmds()
FIU_Fiber = FIU_Fiber_cmds()
modeChange = Mode_Change_cmds()
DM = DM()
sfp = SFP_cmds()
nirc2 = Nirc2_cmds()

nominalFlat = '220405_601'
print('Nominal DM flat for KPIC source:', nominalFlat)

fam_bkgd_pos = [50, 9900]

nolls = {4:'Defocus', 5:'Astig 45', 6:'Astig 0', 7:'Coma Y', 8:'Coma X', 9:'Trefoil Y', 10:'Trefoil X', 11:'Spherical',
         12:'Noll 12', 13:'Noll 13', 14:'Noll 14', 15:'Noll 15', 16:'Noll 16', 17:'Noll 17', 18:'Noll 18', 19:'Noll 19', 20:'Noll 20', 21:'Noll 21', 22:'Noll 22', 23:'Noll 23', 24:'Noll 24'}
# connect to flip mirror keyword
# flip  = ktl.Service("kpic")["PDPONAM"]
# flip.monitor()
serv = ktl.Service("kpic")
serv["PDPONAM"].monitor()
#fib = FEU_Fiber_cmds()
# IS = Telescope_Simulator_cmds()

### These are some parameters that probably don't need to change often
N_order = 7
N_fiber = 1
nskip = 4 ### Only fit every 4 points for speed when using KPIC DRP, still give good enough trace after smoothing
itime = 1500   # for spec
coadds = 1

SAM_off = [3, 4.2]
SAM_on  = [3.292, 4.375]

### For analyzing science data
# gain = 3.03 # e-/ADU
# wv_soln_dat = ascii.read("hip81497_200608_sf2_wvsoln.txt")
# wv_solns = np.array([[ele for ele in row] for row in wv_soln_dat[::-1]])
# k_filt = ascii.read("2massK.txt", names=['wv', 'trans'])
# order_wvs = np.zeros([N_order, 2048])
# for i in range(N_order):
#     pfit  = np.poly1d(wv_solns[i])
#     ### TODO - check that this shape corresponds to the actual orders correctly
#     order_wvs[i] = pfit(np.arange(2048))


def _get_cred2_flux(plot=False, nb_im=50):
    # return 1
    A = np.zeros(nb_im)
    for i in range(A.size):
        ### Get the processed cred2 image
        procim = track.Img.get_data(reform=True)
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

# set optimal reads for NSPEC
def set_nspec_reads(tint):
    if 1400 <= tint < 3000:
        n_reads = 1
    elif 3000 <= tint < 4500:
        n_reads = 2
    elif 4500 <= tint < 6000:
        n_reads = 3
    elif 6000 <= tint <= 60000:
        n_reads = 4
    else:
        n_reads = 16
    # CDS is 2
    if n_reads == 1 and not ktl.read('nspec','sampmode') == 'CDS':
        ktl.write('nspec','sampmode', '2')        
    # MCDS is 2
    elif n_reads > 1 and not ktl.read('nspec','sampmode') == 'MCDS':
        ktl.write('nspec','sampmode', '3')

    ktl.write('nspec','numreads', str(n_reads))
    time.sleep(0.1)

def _acquire_fiber(fibr, verbose=True, toggle_track=False):
    '''
    Change to a new science fiber. If not given a valid fiber number, does nothing
    '''
    if fibr in [1,2,3,4]:
        track.set_goal("scf{}".format(fibr))
        # wait for tracking script to update goal
        while track.get_psf_goal()[3] != fibr:
            time.sleep(.1)

        start_t = time.time()
        # wait for a valid psf to be within .5 pixels in x and y of target
        while True:
            # calculate the dots to print based on how long we've been running
            dots = "."*((int(time.time()) - int(start_t))%3+1) 
            # clear line
            if verbose:
                sys.stdout.write(u"\u001b[2K")
                sys.stdout.write("\rfiber {} | acquiring".format(fibr)+dots)
                sys.stdout.flush()

                # toggle tracking for NCPAs, if stuck here for 5 sec or more
                if toggle_track:
                    cur_t = time.time()
                    if (cur_t-start_t) > 5:
                        track.start_tracking()
                        time.sleep(0.02)
                        track.stop_tracking()

            # get location of psf
            psfx, psfy = track.get_psf_cent()
            # get target location
            goal = track.get_psf_goal()

            # if we're close enough to goal, break
            if track.get_psf_parameters()[0] and abs(psfx - goal[1]) < .5 and abs(psfy - goal[2]) < .5:
                return 1
    else:
        return 0

def _rectify(cube, badpixmap, N_order, N_fiber):
    '''
    Flattens and rectifies a cube for peak finding
    '''
    ### Flatten and get rid of bad pixels
    deg = 4
    ### Badpixmap caussing issues, skipping for now
    # print('cube shape', cube.shape)
    image = np.nansum(cube, axis=0) #*badpixmap
    # plt.imshow(image, norm=LogNorm(vmin=1,vmax=100))
    # plt.show()
    cind = int(image.shape[1]/2)
    nrows = image.shape[0]
    rowinds = np.arange(nrows)
    lind = cind-2
    rind = cind+2
    rectIm = np.empty(np.shape(image)) 
    rectIm[:,cind] = image[:,cind]
    # rectIm[:,cind] = np.nansum(image[:,cind-4:cind+5],axis=1) ### start from the middle, avoid edge effects
    rectCoeffs = np.empty((deg+1,image.shape[1]))
    init = np.zeros(deg+1)
    init[-2] = 1.
    rectCoeffs[:,cind] = init ### coeffs of the identity polynomial

    ### Find full peaks, then retain only the best N_orders*N_fibers 
    if N_fiber==1:
        dist = 120
    else:
        dist = 15
    fullpeaks = signal.find_peaks(np.nansum(image[:,cind-4:cind+5],axis=1), prominence=1e2, distance=dist)
    peaks = fullpeaks[0]

    # plt.plot(rectIm[:,cind])
    # for peak in peaks:
    #   plt.axvline(peak)
    # plt.show()

    proms = fullpeaks[1]['prominences']
    while peaks.size > N_order*N_fiber:
        ### Toss partial traces at detector edge
        if min(peaks) < 40:
            worstpeak = np.argmin(peaks)
        else:
            worstpeak = np.argmin(proms)
        proms = np.delete(proms, worstpeak)
        peaks = np.delete(peaks, worstpeak)
    cpeaks = peaks
    ordcenters = [np.mean(cpeaks[i:i+N_fiber]) for i in range(0,N_order*N_fiber,N_fiber)]

    ### Loop over full image, in chunks
    rectfails = []
    fitinds = []
    while rind < image.shape[1] and lind > 0:   
        ### Get left slice
        lcol = np.nanmedian(image[:,lind-4:lind+5], axis=1)
        lpeaks = signal.find_peaks(lcol, prominence=1e2, distance=dist)
        ### If missing peaks, skip
        if len(lpeaks[0]) < len(cpeaks):
            rectfails.append(lind)
            rectCoeffs[:,lind] = np.nan*np.ones(deg+1)
            rectIm[:,lind] = np.nan*np.ones(rowinds.size)
        else:
            peaks = lpeaks[0]
            ### If we get extra peaks, drop the weakest ones
            if len(peaks) > len(cpeaks):
                proms = lpeaks[1]['prominences']
                while len(peaks) > len(cpeaks):
                    ### Keep the N_fibers peaks closest to the center of each order, throw out rest
                    goodInds = np.empty(N_order*N_fiber)
                    for i in range(len(ordcenters)):
                        center = ordcenters[i]
                        inds = np.argsort([np.abs(peak-center) for peak in peaks])[:N_fiber]
                        goodInds[i*N_fiber:(i+1)*N_fiber] = inds
                    goodInds = np.sort(goodInds).astype(int)
                    peaks = peaks[goodInds]
                lpeaks = peaks
            else: 
                lpeaks = lpeaks[0]
            ### Find polynomial to line up with reference peaks
            # print(lpeaks, cpeaks)
            lcoeffs = np.polyfit(lpeaks, cpeaks, deg=deg)
            rectCoeffs[:,lind] = lcoeffs
            lpoly = np.poly1d(lcoeffs)
            rectIm[:,lind] = np.interp(rowinds, lpoly(rowinds), lcol)
            
        ### Save polynomial coefficients
        fitinds.append(lind)
        lind-=4

        ### Get right slice
        rcol = np.nanmedian(image[:,rind-4:rind+5], axis=1)
        rpeaks = signal.find_peaks(rcol, prominence=1e2, distance=dist)
        ### If missing peaks, skip
        if len(rpeaks[0]) < len(cpeaks):
            rectfails.append(rind)
            rectCoeffs[:,rind] = np.nan*np.ones(deg+1)
            rectIm[:,rind] = np.nan*np.ones(rowinds.size)
        else:
            peaks = rpeaks[0]
            ### If extra peaks, drop weakest ones
            if len(peaks) > len(cpeaks):
                proms = rpeaks[1]['prominences']
                while len(peaks) > len(cpeaks):
                    ### Keep the N_fibers peaks closest to the center of each order, throw out rest
                    goodInds = np.empty(N_order*N_fiber)
                    for i in range(len(ordcenters)):
                        center = ordcenters[i]
                        inds = np.argsort([np.abs(peak-center) for peak in peaks])[:N_fiber]
                        goodInds[i*N_fiber:(i+1)*N_fiber] = inds
                    goodInds = np.sort(goodInds).astype(int)
                    peaks = peaks[goodInds]
                rpeaks = peaks
            else: 
                rpeaks = rpeaks[0]
            ### Find polynomial to line uup with reference peaks
            rcoeffs = np.polyfit(rpeaks, cpeaks, deg=deg)
            rectCoeffs[:,rind] = rcoeffs
            rpoly = np.poly1d(rcoeffs)
            rectIm[:,rind] = np.interp(rowinds, rpoly(rowinds), rcol)
        ### Save polyomial coefficients
        fitinds.append(rind)
        rind+=4

    if len(rectfails) > 0:
        print('Rectification failed in', len(rectfails), 'columns:', np.sort(rectfails))

    rectIm[np.isnan(rectIm)] = 0.

    # plt.imshow(rectIm, norm=LogNorm(vmin=1,vmax=100))
    # plt.show()

    return rectIm, rectCoeffs, cpeaks

def _find_trace(rImage, rCoeffs, goodPeaks, N_order, N_fiber):
    '''
    Given a rectified image, rectification coefficients, peaks, and orders/fibers, finds the traces
    '''
    peaks = goodPeaks
    traces = np.empty((N_order*N_fiber, rCoeffs.shape[1]))
    for i in range(rCoeffs.shape[1]):
        ### Avoid weird things at the edges and bad pixels
        if np.any(np.isnan(rCoeffs[:,i])) or i < 20 or i > 2028:
            traces[:,i] = np.nan*np.ones(N_order*N_fiber)
            continue
        ### Undo the rectification
        for j in range(N_order*N_fiber):
            pcoeff = rCoeffs[:,i]
            pcoeff[np.abs(pcoeff)<1e-16] = 0.
            ### Solve the polynomial 
            poly = np.poly1d(pcoeff)
            try:
                roots = (poly-peaks[j]).roots
            except:
                traces[j,i] = np.nan
                continue
            roots = np.real(roots[np.isreal(roots)])
            root = roots[roots>30] ### Reject bad answers near edges
            if root.size==1: traces[j,i] = float(root)
            elif root.size>1: traces[j,i] = np.min(root-peaks[j])
            else: traces[j,i] = np.nan
        ### Somehow these slip through
        if np.any(np.abs(traces[:,i]-peaks) > 100):
            traces[:,i] = np.nan*np.ones(N_order*N_fiber)
        ### Every now and fibers can switch
        for k in range(len(traces[:,i])-1):
            if traces[k,i] > traces[k+1,i]:
                traces[:,i] = np.nan*np.ones(N_order*N_fiber)
                continue
    ### Get rid of clear outliers that accumulate around edges
    traces[np.abs(traces)>2008] = np.nan
    traces[traces<40] = np.nan

    ### Polynomial nterpolate over NaNs
    for i in range(N_order*N_fiber):
        order = traces[i]
        inds = np.arange(order.size)
        ordfit = np.poly1d(np.polyfit(inds[~np.isnan(order)], order[~np.isnan(order)], deg=4))
        traces[i] = ordfit(inds)

    return traces

def _extract_fluxes(datacube, bkgd_noise, badpixmap, trace_locs, trace_sigs):
    pool = mp.Pool(4)
    flux_lst = []
    for i in range(len(datacube)):
        nframes = str(len(datacube))
        txt = '%3.0f / '%(i+1) + nframes
        sys.stdout.write('\rFlux extraction, frame'+txt)
        sys.stdout.flush()
        image = datacube[i]*badpixmap
        ### Get the diwth and location for each frame
        trace_sig = trace_sigs[i,:,:,:]
        trace_loc = trace_locs[i,:,:,:]

        ### This part was needed to get the flux extraction to work, it breaks if there's only one fiber
        ### illuminated (at least I think that's what's happening). Creates a dummy second fiber.
        if trace_loc.shape[0] == 1:
            arr = 20+trace_loc
            trace_loc = np.append(trace_loc, arr,axis=0)
            trace_sig = np.append(trace_sig, trace_sig, axis=0)
        ### Background traces for DRP
        trace_loc_slit, trace_loc_dark = trace.get_background_traces(trace_loc)
        ### More DRP setup
        trace_loc_wbkg = np.concatenate([trace_loc, trace_loc_slit, trace_loc_dark], axis=0)
        trace_sigmas_wbkg = np.concatenate([trace_sig, trace_sig, trace_sig], axis=0)
        trace_flags = np.array([0, ] * trace_loc.shape[0] + [1, ] * trace_loc.shape[0] + [2, ] * trace_loc.shape[0])
        ### DRP frame processing
        fluxes, errors_extraction, errors_bkgd_only = extraction.extract_flux(image, trace_loc_wbkg, trace_sigmas_wbkg,
                                                                          output_filename=None,
                                                                          img_hdr=None,
                                                                          img_noise=bkgd_noise, fit_background=True,
                                                                          trace_flags=trace_flags, pool=pool,
                                                                          bad_pixel_fraction=0.01,
                                                                          box=True)
        flux_lst.append(fluxes)
    print('')
    int_fluxes = np.nansum(flux_lst, axis=-1)
    ### If averaging < 100 counts/pixel, more likely to run into trace fitting issues
    if np.nanmax(int_fluxes) < 100*N_order*2048:
        print('Warining: Low counts on spec, trace fitting may be unreliable. Check source intensity!')
    return np.asarray(flux_lst)

def _extract_fluxes_hybrid(datacube, bkgd_noise, badpixmap, trace_loc):
    '''
    Does DRP flux extraction using trace locations from the faster trace finding code
    '''
    trace_sigs = np.empty((len(datacube), N_fiber, N_order, 2048))
    if datacube.ndim == 3:
        image = np.nansum(datacube, axis=0)
    else: image = datacube

    image[np.isnan(image)] = 0.0
    image[np.isnan(image)] = 0.0
    _, _, trace_sig, _ =_fit_traces_fast(image, trace_loc, N_order, N_fiber)

    flux_lst = []
    pool = mp.Pool(4)
    for i, image in enumerate(datacube):
        nframes = str(len(datacube))
        txt = '%3.0f / '%(i+1) + nframes
        sys.stdout.write('\rExtracting frame'+txt)
        sys.stdout.flush()

        if trace_loc.shape[0] == 1:
            arr = 20+trace_loc
            trace_loc = np.append(trace_loc, arr,axis=0)
            trace_sig = np.append(trace_sig, trace_sig, axis=0)
        ### Background traces for DRP
        trace_loc_slit, trace_loc_dark = trace.get_background_traces(trace_loc)
        ### More DRP setup
        trace_loc_wbkg = np.concatenate([trace_loc, trace_loc_slit, trace_loc_dark], axis=0)
        trace_sigmas_wbkg = np.concatenate([trace_sig, trace_sig, trace_sig], axis=0)
        trace_flags = np.array([0, ] * trace_loc.shape[0] + [1, ] * trace_loc.shape[0] + [2, ] * trace_loc.shape[0])
        ### DRP frame processing
        fluxes, errors_extraction, errors_bkgd_only = extraction.extract_flux(image, trace_loc_wbkg, trace_sigmas_wbkg,
                                                                          output_filename=None,
                                                                          img_hdr=None,
                                                                          img_noise=bkgd_noise, fit_background=True,
                                                                          trace_flags=trace_flags, pool=pool,
                                                                          bad_pixel_fraction=0.01,
                                                                          box=True)
        flux_lst.append(fluxes)
    print('')
    return np.asarray(flux_lst)

def _extract_fluxes_fast(datacube, trace_loc, plot=False):
    trace_sigs = np.empty((len(datacube), N_fiber, N_order, 2048))
    flux_ests = np.empty((len(datacube), N_fiber, N_order, 2048))
    tot_cnts = np.empty((len(datacube), N_fiber, N_order, 2048))
    trace_locs = np.empty((len(datacube), N_fiber, N_order, 2048))

    if plot:
        collapsedim = np.nansum(datacube,axis=0)
        collapsedim[np.isnan(collapsedim)] = 0.
        plt.title('Stacked cube, with traces')
        plt.imshow(collapsedim, norm=LogNorm(vmin=1,vmax=100))
        for ordloc in trace_loc[0]:
          plt.plot(ordloc, color='r', linestyle='--')
        plt.show()

    for i, image in enumerate(datacube):
        nframes = str(len(datacube))
        txt = '%3.0f / '%(i+1) + nframes
        sys.stdout.write('\rExtracting frame'+txt)
        sys.stdout.flush()

        image[np.isnan(image)] = 0.0
        trace_amp, _, trace_sig, tot_cnt =_fit_traces_fast(image, trace_loc, N_order, N_fiber)
        # trace_locs[i] = trace_loc
        flux_ests[i] = trace_sig*trace_amp*np.sqrt(2*np.pi)
        flux_ests[flux_ests>3e5] = 0.0 ### Reject hot pixel values
        flux_ests[flux_ests<0] = 0.0 ### Reject bad pixel values
        trace_sigs[i] = trace_sig
        tot_cnts[i] = tot_cnt
    print('')
    return flux_ests, tot_cnts

def _fit_traces(datacube, badpixcube, fibers, nskip=4):
    '''
    Use the KPIC DRP to fit each trace for each image in cube. Requires 
    a good guess for fiber positions. Nskip determines how many pixels to skip
    in the fitting, setting it higher speeds things up.
    '''
    # print('badpixshape', badpixcube.shape)
    # plt.imshow(badpixcube[0])
    # plt.show()
    nthreads = 4
    ### Get number of orders from fibers
    norder = len(fibers[list(fibers.keys())[0]])
    _, nx, ny = badpixcube.shape
    ### Fit each trace in the cube
    trace_list = []
    for i in range(len(datacube)):
        nframes = str(len(datacube))
        txt = '%3.0f / '%(i+1) + nframes
        sys.stdout.write('\rTrace fitting, frame'+txt)
        sys.stdout.flush()
        ### Setup for DRP tools - expecting list of images
        image = np.expand_dims(datacube[i],axis=0)
        fiber_list = trace.guess_star_fiber(image, fibers)

        ### Fit traces and smooth using KPIC DRP tools
        trace_calib, residuals = trace.fit_trace(fibers,fiber_list,image[:,:,::nskip],badpixcube[:,:,::nskip],ny,int(nx/nskip),norder,numthreads=nthreads, fitbackground=False)
        # polyfit_trace_calib, smooth_trace_calib = trace.smooth(trace_calib)
        
        ### Save outputs to arrays for plotting
        trace_list.append(trace_calib)
    print('')
    ### Get widths and locations
    trace_list = np.asarray(trace_list)
    # print(trace_list.shape)
    sigs = trace_list[:,:,:,:,1]
    # print(sigs.shape)
    # print(sigs[0,0,1])
    for i in range(len(sigs)):
        try:
            sigs[i] = _fourier_smooth(sigs[i], trace_list.shape[2], trace_list.shape[1])
        except:
            sigs[i] = np.nan*np.ones(np.shape(sigs[i]))
            print('Bad frame', i)

    ### Resample onto full grid if needed
    if sigs.shape[3] != 2048:
        oldshape = sigs.shape
        oldlen = sigs.shape[3]
        oldgrid = np.arange(oldlen)
        newshape = (oldshape[0], oldshape[1], oldshape[2], 2048)
        new_trace_sig = np.empty(newshape)
        newgrid = np.linspace(0, oldlen, num=2048)
        for i in range(newshape[0]):
            for j in range(newshape[1]):
                for k in range(newshape[2]):
                    new_trace_sig[i,j,k] = np.interp(newgrid, oldgrid, sigs[i,j,k])
        sigs = new_trace_sig

    return sigs

def _fit_traces_fast(goodim, traces, N_order, N_fiber):
    rawsig = np.empty((N_fiber, N_order, 2048))
    rawamp = np.empty((N_fiber, N_order, 2048))
    rawtot = np.empty((N_fiber, N_order, 2048))

    if traces.ndim == 2:
        traces = np.expand_dims(traces,axis=0)

    for i in range(N_fiber):
        for j in range(N_order):
            imslice = np.empty((9,2048))
            amps = np.empty(2048)
            for k in range(2048):
                tcen = int(np.round(traces[i,j,k],0))
                cut = np.nansum(goodim[tcen-4:tcen+5,k:k+1], axis=1)

                ### Avoid wackiness at detector edge
                if cut.size != 9 or np.abs(k-1024) > 1020: 
                    imslice[:,k] = np.nan*np.ones(9)
                    amps[k] = np.nan
                else:
                    # bkg = np.nanmedian(np.concatenate((cut[:5], cut[:-5])))
                    imslice[:,k] = cut #-bkg
                    ### Deal with sub-pixel phase by fitting peak
                    p1 = np.polyfit(np.arange(3,6), cut[3:6], 2)
                    pfit = np.poly1d(p1)
                    fitprof = pfit(np.arange(3,6.1,0.1))
                    amps[k] = np.nanmax(fitprof)
                    # plt.plot(cut)
                    # plt.plot(np.arange(3,6.1,0.1), fitprof)
                    # plt.show()
            
            ### Because it's a sampled gaussian, amps will be underestimate on average
            ### Leads to sigs being slightly larger - underestimate resolution
            pshifts = np.abs(np.argmax(imslice, axis=0)-4)
            ### Need to figure out Lorentzian broadening here, thinking it's ~5 percent
            tots = np.nansum(imslice, axis=0)
            sigs = tots/(amps*np.sqrt(2*np.pi))
            amps[pshifts>1.5] = np.nan 
            sigs[pshifts>1.5] = np.nan 
            sigs[sigs>3] = np.nan
            rawamp[i,j] = amps
            rawsig[i,j] = sigs
            rawtot[i,j] = tots
    # print(rawsig)
    ### Fourier smoothing - drop all the high-frequency stuff
    # smoothamp = _fourier_smooth(rawamp, N_order, N_fiber)
    smoothsig = _fourier_smooth(rawsig, N_order, N_fiber)

    return rawamp, rawsig, smoothsig, rawtot

def _fourier_smooth(arr, N_order, N_fiber, cutoff=0.004):
    '''
    Fourier-based smoothing for widths
    '''
    arr[arr<0.5] = np.nan
    smootharr = np.empty(arr.shape)
    try:
        for i in range(N_fiber):
            for j in range(N_order):
                inds = np.arange(arr[i,j].size)
                vals =  arr[i,j]
                badinds = np.isnan(vals) 
                goodinds = inds[~badinds]
                goodvals = vals[~badinds]
                # pdb.set_trace()
                ft = np.fft.fft(goodvals)
                freq = np.fft.fftfreq(goodinds.shape[-1])
                ### This has been tuned manually based on fringe patterns
                ### TODO - replace with a proper window to reduce ringing
                ft.imag[np.abs(freq)>cutoff] = 0.
                ft.real[np.abs(freq)>cutoff] = 0.
                newvals = np.fft.ifft(ft)
                
                smootharr[i,j] = np.interp(inds, goodinds, newvals)
                # plt.plot(vals)
                # plt.plot(smootharr[i,j])
                # plt.show()
        return smootharr
    except ValueError:
        print('Smoothing failed')
        return arr

def _take_backgrounds(itime, coadds, nframes=5, retake=False, Path=None):
    '''
    Take background frames for scan
    '''
    ### Path info
    if Path is None:
        Path = Organization.get_path('Spec')
        Path = list(Path)
        Path[1]=Path[1].split('.')[0]

    ### Check that backgrounds already exist, if they don't force a retake
    bkg_names = []
    for i in range(nframes):
        bkg_names.append(Path[0]+'bkgd_'+str(i)+'_itime'+str(itime)+'.fits')
        if not os.path.isfile(bkg_names[-1]):
            retake = True

    if retake:
        set_nspec_reads(itime)
        ### Stop tracking and move off-slit for backgrounds
        if track.is_tracking():
                track.stop_tracking()
                # Wait for the Tracking loop to complete its last iteration
                time.sleep(2)
        ### Get initial TTM position
        TTM_CP = FAM.get_pos()
        ### Move TTM and wait for new position
        FAM.set_pos(fam_bkgd_pos, block=True)
        proc.save_dark()
        time.sleep(1)

        ### Take backgrounds
        for i in range(nframes):
            ### Take backgrounds, save for DRP 
            spec.TakeImage(Coadds=coadds, tint=itime, Save=True, filename='bkgd_'+str(i)+'_itime'+str(itime))
        FAM.set_pos([TTM_CP[0], TTM_CP[1]], block=True)

    ### Return to original TTM position and restart tracking
    bkgd, smoothed_thermal_noise, badpixmap = background.make_badpixmap(bkg_names)

    track.start_tracking()

    return bkgd, smoothed_thermal_noise, badpixmap

def _take_pd_backgrounds(nreads=10):
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
    # FAM.set_pos([9, 9], block=True)
    FAM.set_pos(fam_bkgd_pos, block=True)
    proc.save_dark()
    time.sleep(1)

    ### Take backgrounds
    with redPM_cmds() as pd:
        bkgd = pd.read_pd(nreads)*1e4

    ### Return to original FAM position and restart tracking
    FAM.set_pos([FAM_CP[0], FAM_CP[1]], block=True)
    if TrackFlag:
        track.start_tracking()

    return bkgd

def _total_flux_fast(datacube, badpixmap, N_order, N_fiber, plotTrace=False):
    '''
    Get the total flux on spec for each image in the cube. Moved here since it's used a bunch
    '''
    ### Rectify for the fast trace finding
    rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber)
    ### Trace fitting
    trace_locs = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)
    
    ### Fast flux extractiong
    if N_fiber == 1:
        trace_locs = np.expand_dims(trace_locs, axis=0)
    _, flux_lst =_extract_fluxes_fast(datacube, trace_locs, plot=plotTrace)

    ### Get total flux in each frame
    int_flux = np.empty(flux_lst.shape[0])
    for i, frame in enumerate(flux_lst):
        int_flux[i] = np.nansum(np.nansum(frame, axis=2), axis=1)[0]
    ### If averaging < 100 counts/pixel, more likely to run into trace fitting issues
    if np.nanmax(int_flux) < 100*N_order*2048:
        print('Warining: Low counts on spec, trace fitting may be unreliable. Check source intensity!')
    return int_flux

def _reduce_drp_byframe(datacube, badpixmap, smoothed_thermal_noise,  N_order, N_fiber, stack=False):
    '''
    Get trace widths and fluxes for each frame in the datacube, using KPIC DRP methods (slow)
    '''
    ### Really wants a cube for this
    
    # badpixmap = np.ones((2048,2048))
    badpixmap = np.ones((2048,2048))
    badpixcube = np.expand_dims(badpixmap, axis=0)
    ### Rectify frame and find trace locations
    # print('Datacube shape', datacube.shape)
    rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber)
    ### Check rectified image if desired
    # plt.imshow(rectIm, norm=LogNorm(vmin=1, vmax=100))
    # plt.show()
    traces = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)
    trace_locs = np.empty((len(datacube), N_fiber, N_order, 2048))
    for i in range(len(datacube)):
        trace_locs[i] = traces

    ### Trace fitting. Currently using KPIC DRP (slow)
    ### This is for the current KPIC DRP trace fitting process
    tracebounds = []
    for order in traces:
        tracebounds.append([int(np.min(order)-20), int(np.max(order)+20)])
    ### Assign to fibers for the DRP
    fibers = {}
    for i in range(N_fiber):
        fibers[i] = tracebounds[i::N_fiber]
    ### Fit trace widths with KPIC DRP
    if not stack:
        trace_sigs = _fit_traces(datacube, badpixcube, fibers, nskip=nskip)
    else:
        trace_sigs = []
        trace_sig = _fit_traces(np.expand_dims(np.nansum(datacube,axis=0),axis=0), badpixcube, fibers, nskip=nskip)
        for i in range(len(datacube)):
            trace_sigs.append(trace_sig[0])
        trace_sigs = np.asarray(trace_sigs)
    # ## Do flux extraction with KPIC DRP tools, box option
    flux_lst = _extract_fluxes(datacube, smoothed_thermal_noise, badpixmap, trace_locs, trace_sigs)

    return trace_sigs, trace_locs, flux_lst

def Scam_Reduced(nb_TC_im = 10, Method = 'GaussFit'):
    pass
    ''' -----------------------------------------------------------------------
    ----------------------------------------------------------------------- '''
    # Take image(s) with the Tracking Camera before to take Scam image
    # C_Im_I = np.mean(cam.grab_n(nb_TC_im)[0].data, 0)
        
    # # Take image(s) with the FEU detector.
    # Scam_tint = 656
    # scam.TakeImage(Coadds = 1,tint = Scam_tint)
        
    # # Store the reduced image in the dedicated cube of images
    # S_Im = scam.im[1,:,:]/Scam_tint
        
    # # Take image(s) with the Tracking Camera after to take Scam images
    # C_Im_E = np.mean(cam.grab_n(nb_TC_im)[0].data, 0)

    # # Compute the flux of the PSF on the Cred2 (Before)
    # try: 
    #     PSF_Ini = PSF_Finder(C_Im_I,style = 'GaussFit')   
    # except: 
    #     PSF_Ini = False
                    
    # # Compute the flux of the PSF on the Cred2 (After)
    # try: 
    #     PSF_End = PSF_Finder(C_Im_E,style = 'GaussFit')
    # except: 
    #     PSF_End = False
    # # print(PSF_Ini)
    
    # # Check if both Cred2 flux measurements are valid
    # if (PSF_Ini == False and PSF_End == False): 
    #     PSF_Ini,PSF_End = [1,0,0,0,0,0,0],[1,0,0,0,0,0,0]
    # elif PSF_Ini == False: 
    #     PSF_Ini = PSF_End
    # elif PSF_End == False: 
    #     PSF_End = PSF_Ini

    # # Compute Cred PSF Parameters
    # Cred_PSF    = np.round((PSF_Ini[0]+PSF_End[0])/2.,3)

    # if Method == 'Photometry':
    #     # Compute Scam PSF Parameters 
    #     try: 
    #         Scam_PSF = scam.Photometry(S_Im,3)[0]
    #     except: 
    #         Scam_PSF = np.zeros([3])
    
    # elif Method == 'GaussFit':
    #     # 
    #     try: 
    #         Scam_PSF = PSF_Finder(S_Im,style = 'GaussFit')[0]   
    #     except: 
    #         Scam_PSF = np.zeros([7]) 
    #     if np.all(Scam_PSF) == False: 
    #         Scam_PSF = np.zeros([7])            
    
    # return S_Im, C_Im_I, C_Im_E, Scam_PSF, Cred_PSF, Scam_tint

def PSM_Scan_XY(xstart=None, xstop=None, xstep=2., ystart=None, ystop=None, ystep=2., Scam_tint=1500, fiber=2, usavename='', plot=True, keepTracking=False):
    ''' -----------------------------------------------------------------------
    Scan the Zabers in x/y
    ----------------------------------------------------------------------- '''

    # Get Zaber initial positions
    PSMX_ini, PSMY_ini = PSM.get_pos()

    if xstart is None:
        xstart = PSMX_ini-xstep*4
    if xstop is None:
        xstop = PSMX_ini+xstep*4
    if ystart is None:
        ystart = PSMY_ini-ystep*4
    if ystop is None:
        ystop = PSMY_ini+ystep*4


    # Set the FEU TTM so the PSF is out of the slit
    SAMX_ini, SAMY_ini = SAM.get_pos()
    # Verify if tracking loop active.
    #Track_Flag = track.is_tracking()
    goal0 = track.get_goal()[0]
    try:
        SAM.set_pos("off_slit")

        # Prepare scanned position
        List_pos_x = np.arange(xstart, xstop+xstep/2., xstep)
        List_pos_y = np.arange(ystart, ystop+ystep/2., ystep)
        # print(List_pos_x, List_pos_y)

        # Compute the total number of steps
        nbsteps = np.size(List_pos_x)*np.size(List_pos_y)

        # Do not turn off tracking - edited 2022, Oct 8
        # if Track_Flag:
        #     track.stop_tracking()
        #     # Wait for the Tracking loop to complet its last iteration
        #     time.sleep(20 * cam.get_tint())
            # print('Tracking loop gain set to 0 for the scan.')
        
        # Get the current position of the TTM
        TTM_CP = FAM.get_pos()

        # Move the TTM to take backgrounds
        track.stop_tracking()
        FAM.set_pos(fam_bkgd_pos, block = True)

        # Take Scam calibration images
        scam.TakeImage(tint=Scam_tint)
        scam_bkgd = scam.im[0]

        # Return to initial TTM position
        FAM.set_pos([TTM_CP[0],TTM_CP[1]], block=True)
        track.start_tracking()
        # go to specified fiber
        _acquire_fiber(fiber)

        # Init the iteration number
        ite_nb  = 0
        # Prepare the cube of image to store raw data (Scam)
        S_Cube = np.zeros([nbsteps,np.shape(scam.im)[1],np.shape(scam.im)[2]])
        
        # Get the initial time
        Time_ini = time.time()

        # Instancy the variable Results.
        # It will contain:
        #    - Flux
        #    - Position of the Zaber X
        #    - Position of the Zaber Y
        Results = np.zeros([3,np.size(List_pos_x),np.size(List_pos_y)])
        Scam_images = np.zeros([np.size(List_pos_x),np.size(List_pos_y),256, 256])
        # cred2fluxes = np.zeros([np.size(List_pos_x),np.size(List_pos_y)])
        for pos_x in List_pos_x:
            # using move parameter means we only move the axis we want
            PSM.set_pos(target=[pos_x, PSMY_ini], block=True, move=[1,0])
            for pos_y in List_pos_y:
                # Move Zabers to position
                PSM.set_pos(target=[pos_x, pos_y], block=True, move=[0,1])
                
                # Compute iteration number in x and y directions
                tmp_x = np.where(List_pos_x == pos_x)[0]
                tmp_y = np.where(List_pos_y == pos_y)[0]
                
                ### Take cred2 data
                cred2flux = _get_cred2_flux()
                # Take images with scam and reduce the data
                scam.TakeImage(tint=Scam_tint)
                Scam_images[tmp_x,tmp_y,:,:] = scam.im[0]-scam_bkgd
                #tmp_im = medfilt(scam.im[0]-scam_bkgd,3)
                #x,y = np.where(tmp_im == np.max(tmp_im))
                #x = int(x[0])
                #y = int(y[0])
                #### Get counts in area around max flux on SCAM
                #flux = np.sum(tmp_im[x-3:x+4,y-3:y+4])
                #Results[0,tmp_x,tmp_y] = flux #/cred2flux
                ## Position of the Zaber X
                #Results[1,tmp_x,tmp_y] = pos_x
                ## Position of the Zaber Y
                #Results[2,tmp_x,tmp_y] = pos_y  

                # Compute average time per iteration
                avr_time = (time.time()-Time_ini)/(ite_nb+1.)
                # Compute time before end of scan
                timeleft = round((nbsteps-ite_nb)*avr_time,0)
                # Increment ite number
                ite_nb += 1
                # Prepare information to print
                text  = 'Remaining Time %05.d sec' %(timeleft)
                text += 'X = %05.3f, Y = %05.3f'%(pos_x, pos_y)
                text += ' -- %03d/%03d ite left' %(nbsteps-ite_nb,nbsteps)
                #text += ' -- Flux = %06.3f' %(Results[0,tmp_x,tmp_y])
                # Print the information in the terminal
                sys.stdout.write('\r PSM_Scan_XY: ' + text)
                sys.stdout.flush()


        tot_img = np.sum(Scam_images, axis=(0,1))
        tmp_im = medfilt(tot_img, 3)
        y0, x0 = np.where(tmp_im == np.max(tmp_im))
        x0 = int(x0[0])
        y0 = int(y0[0])
        for pos_x in List_pos_x:
            for pos_y in List_pos_y:
                # Compute iteration number in x and y directions
                tmp_x = np.where(List_pos_x == pos_x)[0][0]
                tmp_y = np.where(List_pos_y == pos_y)[0][0]
            
                tmp_im = Scam_images[tmp_x,tmp_y,:,:]
                # find peak pixel
                tmp_im_filt = medfilt(tmp_im, 3)
                searchrad = 15
                cutout = tmp_im_filt[y0-searchrad:y0+searchrad+1, x0-searchrad:y0+searchrad+1]
                y, x = np.unravel_index(np.argmax(cutout), cutout.shape)
                y += (y0 - searchrad)
                x += (x0 - searchrad)
                
                ### Do simple circular aperture photometry
                ycoord, xcoord = np.indices(tmp_im.shape)
                dist_from_source = np.sqrt((xcoord - x)**2 + (ycoord - y)**2)
                in_aper = np.where(dist_from_source < 8)
                bkgd_aper = np.where((dist_from_source > 15) & (dist_from_source < 20))
                flux_aper = np.sum(tmp_im[in_aper])
                bkgd_level = np.median(tmp_im[bkgd_aper])
                flux = flux_aper - np.size(in_aper) * bkgd_level
                
                ### Get counts in area around max flux on SCAM
                # flux = np.sum(tmp_im[y-12:y+13,x-12:x+13])
                Results[0,tmp_x,tmp_y] = flux #/cred2flux
                # Position of the Zaber X
                Results[1,tmp_x,tmp_y] = pos_x
                # Position of the Zaber Y
                Results[2,tmp_x,tmp_y] = pos_y  
                print(flux, pos_x, pos_y, tmp_x, tmp_y)
                
                
        # Return to initial TTM  and Zaber positions
        FAM.set_pos([TTM_CP[0],TTM_CP[1]], block=True)
        SAM.set_pos([SAMX_ini, SAMY_ini], block=True)
        PSM.set_pos([PSMX_ini, PSMY_ini], block=True)
        
        #if Track_Flag: track.start_tracking()
        
        # Get path where data has to be saved
        Path, filename = Organization.get_path('Scan_FEU_XY')

        # --- Save XXXX
        # Create the name of the data
        fullname = Path + filename + '_Results.fits'
        # Create a Header Data Unit (HDU) based on the data to save.
        hdu = fits.PrimaryHDU(Results)
        hdu.header['fiber'] = goal0
        # Save the data
        hdu.writeto(fullname, overwrite = True)

        # --- Save XXXX
        # Create the name of the data
        fullname = Path + filename + '_Scam_Images.fits'
        # Create a Header Data Unit (HDU) based on the data to save.
        hdu = fits.PrimaryHDU(Scam_images)
        hdu.header['fiber'] = goal0
        # Save the data
        hdu.writeto(fullname, overwrite = True)

        # Extract the flux
        Flux        = Results[0,:,::-1]
        # Extract the Zaber position in both directions
        Zaber_pos_X = Results[1,:,:]
        Zaber_pos_Y = Results[2,:,:]
       
        ### Plot output - check X/Y!!!
        Flux = (Flux - np.min(Flux))/np.max(Flux - np.min(Flux))
        # 
        Rx, Ry = np.meshgrid(List_pos_x, List_pos_y)
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
        p,u = opt.curve_fit(Gaussian2D, dim, Vec, p0 = SP)

        fitmod = np.reshape(Gaussian2D(dim,p[0],p[1],p[2],p[3],p[4],p[5],p[6]),dim) 
        residual = Flux - fitmod
        #
        offset_x = fx(p[1])
        offset_y = fy(p[2])
        #
        print('Gaussfit maximum: %8.6f' %(p[0]))
        print('Best FEU TTM position (fixed): X = %06.2f -- Y = %06.2f' %(offset_x,offset_y))

        if plot:
        	### Check if axes are flipped. Also if it's not moving right
            fig = _plot_fam_scan(Rx, Ry, Flux.T, fitmod.T, residual.T, offset_x, offset_y, (Path,filename), usavename)
            fig.show()

        return offset_x, offset_y

    except (Exception, KeyboardInterrupt) as e:
        print('Exception encountered, resetting initial state')
        FAM.set_pos([TTM_CP[0],TTM_CP[1]], block=True)
        SAM.set_pos([SAMX_ini, SAMY_ini], block=True)
        PSM.set_pos([PSMX_ini, PSMY_ini], block=True)
        #if Track_Flag: track.start_tracking()
        print('Positions reset')
        # Get path where data has to be saved
        Path, filename = Organization.get_path('Scan_FEU_XY')

        # --- Save XXXX
        # Create the name of the data
        fullname = Path + filename + '_Results.fits'
        # Create a Header Data Unit (HDU) based on the data to save.
        hdu = fits.PrimaryHDU(Results)
        hdu.header['fiber'] = goal0
        # Save the data
        hdu.writeto(fullname, overwrite = True)

        # --- Save XXXX
        # Create the name of the data
        fullname = Path + filename + '_Scam_Images.fits'
        # Create a Header Data Unit (HDU) based on the data to save.
        hdu = fits.PrimaryHDU(Scam_images)
        hdu.header['fiber'] = goal0
        # Save the data
        hdu.writeto(fullname, overwrite = True)
        # print(e)
        raise e
        return np.nan, np.nan 

def _take_focus_scan(start=8.8, stop=9.5, step=0.08, bkgd=None, Path=None, savename='', fibers=None):
    '''
    Function to take focus scan images and save a background-subtracted cube at the specified focus points.
    Will reuse background frames if available, or take new ones. Returns the foci used in the scan 
    as well as a numpy cube of the 2D scan images. savename is a string to add to the end of the filename,
    if desired.
    '''
    F_ini = PSM.get_foc_pos(foc=2)[0]
    Track_Flag = track.is_tracking()
    if Track_Flag:
        goal0 = track.get_goal()[0]
    try:
        if bkgd is None:
            bkgd, _, _, = _take_backgrounds(itime, coadds, nframes=5)
        if Path is None:
            Path = Organization.get_path('Spec')
            Path = list(Path)
            Path[1] =Path[1].split('.')[0]
        if fibers is None:
            fibers = [goal0]

        set_nspec_reads(itime)

        datacubes = []
        for fibr in fibers:
            print('Fiber', fibr)
            # add extra line for \033[F\033[K to delete
            print()
            _ = _acquire_fiber(fibr)
            foci = np.arange(start, stop+step/2., step)
            datacube = np.empty((foci.size, 2048, 2048))
            ipos = PSM.get_foc_pos(foc=2)[0]
            for i in range(foci.size):
                try:
                    ### Move FEU
                    PSM.set_foc_pos(foci[i], foc=2, block=True)

                    ### Get and write current position to terminal
                    tmp = PSM.get_foc_pos(foc=2, update = False)[0]
                    txt = 'focus = %5.2f' %(tmp)
                    print('\033[F\033[KCurrent position: ' + txt)

                    ### Take an image, get data-background
                    spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='focus_scan_step'+str(i))
                    datacube[i] = spec.im[0].T[:,::-1] - bkgd
                except:
                    datacube[i] = np.zeros((2048,2048))
            ### Save datacube 
            usavename = '_fiber'+str(fibr)+savename
            hdu = fits.PrimaryHDU(datacube)
            fstr = '%5.2f | '*foci.size
            hdu.header['Foci'] = fstr%tuple(foci)
            hdu.header['fiber'] = fibr
            hdu.writeto(Path[0] + Path[1] + '_focus_scan_cube'+usavename+'.fits', overwrite=True)
            print('')
            ### Set FEU back to initial position
            try:
                PSM.set_foc_pos(F_ini, foc=2, block=True)
            except:
                print('PSM focus failure, position: ', ipos)
            datacubes.append(datacube)

        if len(datacubes)==1:
            return foci, datacubes[0]
        else:
            return foci, datacubes

    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        PSM.set_foc_pos(F_ini, foc=2, block=True)
        if Track_Flag: 
            track.set_goal(goal0)
            track.start_tracking()
        print('Reset complete')
        # print(e)
        raise e
        return np.asarray([np.nan]), np.nan*np.zeros((2048,2048))

def _analyze_focus_scan(datacube=None, foci=None, badpixmap=np.ones((2048,2048)), smoothed_thermal_noise=None, fibers=None, plot=True, Path=None, savename=''):
    '''
    Makes focus scan plots for the given datacube. If plot is True, will save plots to specified Path.
    If foci not given as an array, will compute as arbitrary steps. Separate function so that the 
    analysis can be easily run on whatever datacube you want, not just when the data is first taken.
    If datacube not given, prompts user for entry.
    '''
    if Path is None:
        Path = ('', '')
    ### Use currently selected fiber if unspecified
    if fibers is None:
        goal0 = track.get_goal()[0]
        fibers = [goal0]    ### Use a list of cubes for multiple fibers
    if datacube is not None:
        datacube = np.asarray(datacube)
    datacubes = []
    if datacube is None:
        for i in range(len(fibers)):
            cname = input('Enter datacube name, fiber '+str(fibers[i])+' >> ').replace('\'', '').replace(' ', '')
            fname = Path[0]+Path[1]+cname
            if not os.path.exists(fname):
                print('Invalid datacube filename!')
                print(fname)
                return -1
            else:
                datacubehdu = fits.open(fname)
                datacubes.append(datacubehdu[0].data)
                foci = np.asarray(datacubehdu[0].header['Foci'].split('|')[:-1],dtype=float)
                Path = ('', fname.split('_focus')[0])

    elif datacube.ndim==3:
        datacubes.append(datacube)
    else:
        datacubes = datacube
    datacubes = np.asarray(datacubes)
    # print('shape', datacubes.shape)
    ### Get other things
    if smoothed_thermal_noise is None:
        smoothed_thermal_noise=np.zeros((2048,2048))
    if foci is None:
        foci = np.arange(len(datacube))
        ustr = 'steps'
    else: ustr = 'mm'

    # foci = foci[::3]

    best_foci = []
    alt_foci = []
    for i, datacube in enumerate(datacubes):
        fibr = fibers[i]
        print('Fiber', fibr)
        usavename = '_fiber'+str(fibr)+savename

        # datacube = datacube[::3]
        # print('pre-removal', datacube.shape)
        goodims = []
        for inum in range(len(foci)):
            # print('frame number', inum)
            # print(np.abs(np.nansum(np.nansum(datacube[inum],axis=0),axis=0)))
            if np.abs(np.nansum(np.nansum(datacube[inum],axis=0),axis=0)) > 1e-14:
                goodims.append(inum)
        print('Good images', goodims)
        goodims = np.asarray(goodims,dtype=int)
        # print(goodims)
        foci = foci[goodims]
        datacube = datacube[goodims,:,:]
        # for im in datacube:
        #   plt.imshow(im, Norm=LogNorm(vmin=1, vmax=100))
        #   plt.show()
        # print('post-removal', datacube.shape)

        trace_sigs, trace_locs, flux_lst = _reduce_drp_byframe(datacube, badpixmap, smoothed_thermal_noise, N_order, N_fiber)

        ### Average width/flux over each order
        widths = np.empty((foci.size, N_order))
        fluxes = np.empty((foci.size, N_order))
        for j in range(len(trace_sigs)):
            # print('cube shape',datacube.shape)
            if np.abs(np.nansum(np.nansum(datacube[j],axis=0),axis=0)) < 1e-14:
                widths[j] = np.nan 
                fluxes[j] = np.nan
                continue
            width = trace_sigs[j,0,:,:]
            width = np.nanmedian(width, axis=1)*2.355 ### To FWHM
            widths[j] = width
            flux = np.nanmedian(flux_lst[j,0], axis=1)
            fluxes[j] = flux
        widths = widths.transpose()
        fluxes = fluxes.transpose() 
        fluxes[np.abs(fluxes)>1e4] = np.nan

        ### Polyfit for each order, interpolate fluxes
        xfit = np.arange(foci[0], foci[-1], 0.01)
        fitwidths = np.empty((N_order, xfit.size))
        fitfluxes = np.empty((N_order, xfit.size))
        for j, order in enumerate(widths):
            pcos = np.polyfit(foci[~np.isnan(order)], order[~np.isnan(order)], 4)
            pfit = np.poly1d(pcos)
            fitwidths[j] = pfit(xfit)
            fitfluxes[j] = np.interp(xfit, foci[~np.isnan(order)], fluxes[j][~np.isnan(order)])

        ### Find options near 2.5 pixel FWHM
        meds = np.median(fitwidths,axis=0)
        diffs = np.abs(meds-2.5)
        bestpts = np.argsort(diffs)
        badpts = np.argwhere(np.min(fitwidths,axis=0)<2.)
        ### Want the two best options
        options = np.empty(2, dtype=np.int32)
        k,q = 0,0
        while k < options.size:
            if bestpts[q] not in badpts:
                if np.abs(bestpts[q] - options[k-1]) > 5:
                    options[k] = int(bestpts[q])
                    k+=1
                q+=1
            else:
                q+=1
        ### Use flux as the tiebreaker
        best_focus = xfit[options[np.argmax(fitfluxes[2][options])]]
        alt_focus = xfit[options[np.argmin(fitfluxes[2][options])]]
        best_foci.append(np.round(best_focus,3))
        alt_foci.append(np.round(alt_focus, 3))
        print('Best focus:', np.round(best_focus,3), ustr)
        print('Alternate focus:', np.round(alt_focus,3), ustr)

        ### Make plots if wanted
        if plot:
            mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=pl.cm.plasma(np.linspace(0,1, N_order)[::-1]))
            for i in range(len(widths)):
                plt.scatter(foci, widths[i], label='Order '+str(i))
            for i in range(len(widths)):
                plt.plot(xfit, fitwidths[i])
            plt.axvline(x=best_focus, c='r', linestyle='--', label='Best Focus')
            plt.axvline(x=alt_focus, c='g', linestyle='--', label='Alt Focus')
            plt.axhline(y=2.0, c='k', linestyle='--')
            plt.axhline(y=3.0, c='k', linestyle=':')
            plt.ylim(bottom=1.5, top=5)
            plt.ylabel('Trace FWHM [pix]')
            plt.xlabel('Focus Step ['+ustr+']')
            if fibr !=-1: plt.title('Fiber '+str(fibr))
            plt.legend()
            plt.savefig(Path[0]+Path[1]+'_width_v_focus'+usavename+'.png', bbox_inches='tight')
            plt.show()

            ### Plot flux in each order vs step in scan
            for i in range(len(fluxes)):
                plt.scatter(foci, fluxes[i], label='Order '+str(i))
            for i in range(len(fluxes)):
                plt.plot(xfit, fitfluxes[i])
            plt.axvline(x=best_focus, c='r', linestyle='--', label='Best Focus: '+str(best_focus))
            plt.axvline(x=alt_focus, c='g', linestyle='--', label='Alt Focus: '+str(alt_focus))
            plt.ylabel('Pixels')
            plt.xlabel('Focus ['+ustr+']')
            if fibr !=-1: plt.title('Fiber '+str(fibr))
            plt.legend()
            plt.savefig(Path[0]+Path[1]+'_flux_v_focus'+usavename+'.png', bbox_inches='tight')
            plt.show()
        ### Save stuff for later plotting
        np.savez(Path[0]+Path[1]+'_focus_scan_data'+usavename+'.npz', foci=foci, widths=widths, fluxes=fluxes, locs=trace_locs, rawsigs=trace_sigs, rawfluxes=flux_lst)

    if len(best_foci) == 1:
        return best_foci[0], alt_foci[0]
    else:
        return best_foci, alt_foci

def focus_scan(start=8.8, stop=9.5, step=0.08, plot=True, fibers=None):
    '''
    Scan between given focus values for the FEU position, inclusive bounds
    '''
    # ### Get paths
    Path = Organization.get_path('Spec')
    Path = list(Path)
    Path[1] =Path[1].split('.')[0]
    ### Take backgrounds
    ### For some reason the badpixmap is bad on the Keck machine, but fine locally
    bkgd, smoothed_thermal_noise, badpixmap = _take_backgrounds(itime, coadds, nframes=5)
    badpixmap = np.ones((2048,2048))

    ### Take scan data
    foci, datacube = _take_focus_scan(start=start, stop=stop, step=step, bkgd=bkgd, Path=Path, fibers=fibers)

    ### Analyze the scan 
    best_focus, alt_focus = _analyze_focus_scan(datacube=datacube, foci=foci, fibers=fibers, badpixmap=badpixmap, smoothed_thermal_noise=smoothed_thermal_noise, plot=plot, Path=Path)

    return best_focus

def _take_fam_spec_scan(start=-100, stop=100, nsteps=20, start_pos=None, bkgd=None, Path=None, savename='', save=True):
    '''
    Perform a 2D FAM scan with NIRSPEC as the detector (instead of PD)
    NOTE: start and stop are in FAM units. nsteps is number of steps in the scan in both axes. start_pos should be a 2 element array or list with X,Y coordinates for the center of the scan or None to use the current FAM coords when the function is called.

    Saves a background-subtraced cube to path (NFIU Spec path if none), with optional 
    additional string savename. 
    '''
    pos0 = FAM.get_pos()
    Track_Flag = track.is_tracking()
    try:
        if bkgd is None:
            bkgd, _, _, = _take_backgrounds(itime, coadds, nframes=5)
        if Path is None:
            Path = Organization.get_path('Spec')
            Path = list(Path)
            Path[1] =Path[1].split('.')[0]
        if start_pos is None:
            start_pos = pos0
        if Track_Flag:
            # Turn off tracking to work in FAM coords
            track.stop_tracking()
        
        fam_dels = np.linspace(start, stop, nsteps)
        fam_delxs = fam_dels+start_pos[0]
        fam_delys = fam_dels+start_pos[1]

        datacubes = []
        set_nspec_reads(itime)

        print('Center Position: {}'.format(start_pos))
        datacube = np.empty((nsteps, nsteps, 2048, 2048))
        # sample frame for size
        im = cam.grab_n(1)[0].data
        crd2cube = np.empty((nsteps, nsteps, im.shape[1], im.shape[2]))
        ### Scan over pts
        for xind, fam_delx in enumerate(fam_delxs):
            for yind, fam_dely in enumerate(fam_delys):
                newpos = np.array([fam_delx, fam_dely])

                ### Move FAM
                FAM.set_pos(newpos, block=True)
                tmp = FAM.get_pos(update = False)
                txt = 'X = %5.2f - Y = %5.2f' %(tmp[0], tmp[1])
                sys.stdout.write('\rFAM Position: ' + txt)
                sys.stdout.flush()

                ### Take CRED2 image
                im = cam.grab_n(1)[0].data
                ### Take an image, get data-background
                spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='slit_scan_step'+str(i))
                datacube[xind, yind] = spec.im[0].T[:,::-1] - bkgd
                crd2cube[xind, yind] = im
        
        if save:
            # Save Spec data
            hdu = fits.PrimaryHDU(datacube)
            hdu.header['start'] = start
            hdu.header['stop'] = stop
            hdu.header['nsteps'] = nsteps
            hdu.header['strtposX'] = start_pos[0]
            hdu.header['strtposY'] = start_pos[1]
            tstamp = str(time.time())
            hdu.writeto(Path[0] + Path[1] + '_fam_scan_cube_'+tstamp+'.fits', overwrite=True)
            # Save CRED2 data
            hdu = fits.PrimaryHDU(crd2cube, header=cam.grab_n(1)[0].header)
            hdu.header['start'] = start
            hdu.header['stop'] = stop
            hdu.header['nsteps'] = nsteps
            hdu.header['strtposX'] = start_pos[0]
            hdu.header['strtposY'] = start_pos[1]
            hdu.writeto(Path[0] + Path[1] + '_fam_scan_cube_cred2_'+tstamp+'.fits', overwrite=True)

        print('')
        ## Set FAM back to initial position
        if Track_Flag:
            # Start tracking again
            track.start_tracking()
        else:
            # Set to first position
            FAM.set_pos(pos0, block=True)
        datacubes.append(datacube)

        ### Return either the datacube, or the list of datacubes 
        if len(datacubes) == 1:
            return fam_delx, fam_dely, datacubes[0], crd2cube
        else:
            return xfam_delx, fam_dely, datacubes, crd2cube
    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        FAM.set_pos(pos0, block=True)
        if Track_Flag: 
            track.start_tracking()
        else:
            FAM.set_pos(pos0, block=True)
        print('Reset complete')
        raise e
        assert 1==0
        return np.asarray([np.nan]), np.nan*np.zeros((2048,2048))


def _take_slit_scan(start=3.275, stop=3.325, step=0.005, bkgd=None, Path=None, savename='', fibers=None, save=True,axis='x'):
    '''
    Steps through FEU TTM positions between start and stop (inclusive) with step size step.
    Saves a background-subtraced cube to path (NFIU Spec path if none), with optional 
    additional string savename. Will do for each fiber in the list of fibers given, or the
    current fiber if None
    '''
    pos0 = SAM.get_pos()
    
    Track_Flag = track.is_tracking()
    if Track_Flag:
        goal0 = track.get_goal()[0]
    else: goal0 = 0
    try:
        if bkgd is None:
            bkgd, _, _, = _take_backgrounds(itime, coadds, nframes=5)
        if Path is None:
            Path = Organization.get_path('Spec')
            Path = list(Path)
            Path[1] =Path[1].split('.')[0]
        if fibers is None:
            fibers = [goal0]

        set_nspec_reads(itime)

        datacubes = []
        for fibr in fibers:
            print('Fiber', fibr)
            _ = _acquire_fiber(fibr)

            xpts = np.arange(start, stop+step/2., step)
            datacube = np.empty((xpts.size, 2048, 2048))
            ipos = SAM.get_pos()
            ### Scan over pts
            for i in range(xpts.size):
                ### Move FEU TTM
                if axis=='x':
                    SAM.set_pos([xpts[i], ipos[1]], block=True, move=[1,0])
                elif axis=='y':
                    SAM.set_pos([ipos[0], xpts[i]], block=True, move=[0,1])
                else:
                    raise ValueError('Invalid axis command!')

                tmp = SAM.get_pos(update = False)
                txt = 'X = %5.4f - Y = %5.4f' %(tmp[0], tmp[1])
                sys.stdout.write('\rFEU TTM Position: ' + txt)
                sys.stdout.flush()

                ### Take an image, get data-background
                spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='slit_scan_step'+str(i))
                datacube[i] = spec.im[0].T[:,::-1] - bkgd
            
            if save:
                usavename = '_fiber'+str(fibr)+savename
                hdu = fits.PrimaryHDU(datacube)
                xstr = '%5.4f | '*xpts.size
                hdu.header['xpts'] = xstr%tuple(xpts)
                hdu.header['fiber'] = fibr
                hdu.writeto(Path[0] + Path[1] + '_slit_scan_cube'+usavename+'.fits', overwrite=True)
            print('')
            
            ## In the end, put FEU back to first scan location - IMPORTANT for Repeatability
            if axis == 'y':
                # first set FEU back to initial position
                SAM.set_pos(ipos, block=True, move=[0,1])
                # set to the first scan point
                SAM.set_pos([ipos[0], xpts[0]], block=True, move=[0,1])
                
            elif axis == 'x':
                SAM.set_pos(ipos, block=True, move=[1,0])
                SAM.set_pos([xpts[0], ipos[1]], block=True, move=[1,0])
            
            datacubes.append(datacube)

        ### Return either the datacube, or the list of datacubes for each fiber
        if len(datacubes) == 1:
            return xpts, datacubes[0]
        else:
            return xpts, datacubes
    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        SAM.set_pos(pos0, block=True)
        if Track_Flag: 
            track.set_goal(goal0)
            track.start_tracking()
        print('Reset complete')
        # print(e)
        raise e
        assert 1==0
        return np.asarray([np.nan]), np.nan*np.zeros((2048,2048))

def _analyze_slit_scan(datacube=None, xpts=None, badpixmap=np.ones((2048,2048)), plot=True, Path=None, fibers=None, savename='', plotTrace=False, N_order_local=7):
    '''
    Gets total flux and makes plot for the specified slit scan. If datacube is none, prompts user for file
    '''
    if datacube is not None:
        datacube = np.asarray(datacube)
    ### Getting inputs
    if Path is None:
        Path = ('', '')
    if fibers is None:
        goal0 = track.get_goal()[0]
        fibers = [goal0]
    ### Use a list of cubes for multiple fibers
    datacubes = []
    if datacube is None:
        for i in range(len(fibers)):
            cname = input('Enter datacube name, fiber '+str(fibers[i])+' >> ').replace('\'', '').replace(' ', '')
            fname = Path[0]+Path[1]+cname
            if not os.path.exists(fname):
                print('Invalid datacube filename!')
                print(fname)
                return -1
            else:
                datacubehdu = fits.open(fname)
                datacubes.append(datacubehdu[0].data)
                xpts = np.asarray(datacubehdu[0].header['xpts'].split('|')[:-1],dtype=float)
                Path = ('', fname.split('_slit')[0])
    elif datacube.ndim==3:
        datacubes.append(datacube)
    else:
        datacubes = datacube
    datacubes = np.asarray(datacubes)
    ### Xpts are steps if not given
    if xpts is None:
        xpts = np.arange(len(datacube))
        ustr = 'step'
    else: ustr = 'pos'

    transmissions = np.ones(xpts.size)
    for i, datacube in enumerate(datacubes):
        fibr = fibers[i]
        print('Fiber '+str(fibr))
        usavename = '_fiber'+str(fibr)+savename
        ### Estimate total flux
        int_flux = _total_flux_fast(datacube, badpixmap, N_order_local, N_fiber, plotTrace=plotTrace)
        transmissions*=int_flux/np.nanmax(int_flux)

        ### Output and plot
        best_pos = np.round(xpts[np.argmax(int_flux)],4)
        print('Best FEU TTM position:', best_pos, ustr)
        np.savez(Path[0]+Path[1]+'_slit_scan_data'+usavename+'.npz', int_flux=int_flux)
        print('It is recommended to take a few more spec frames after moving to new SAM position')
        if plot:
                plt.plot(xpts, int_flux, label='Fiber '+str(fibr))  
    if plot:
        plt.axvline(x=xpts[np.argmax(transmissions)], color='r', linestyle='--',label='Best FEU TTM Position')
        plt.text(xpts[np.argmax(transmissions)]+0.005, 0.75*np.max(int_flux), str(np.round(xpts[np.argmax(transmissions)],4)))
        plt.ylabel('Total Flux on Spec')
        plt.xlabel('FEU TTM x position ['+ustr+']')
        plt.legend()
        if len(fibers) > 1:
            usavename = '_fibers'
            for fibr in fibers:
                usavename = usavename+str(fibr)
            usavename = usavename+savename
        plt.savefig(Path[0]+Path[1]+'_slit_scan'+usavename+'.png', bbox_inches='tight')
        plt.show()
        
    ## Whether to move
    # confirm = input('Update position of SAM to '+str(best_pos)+ '? [Y/n] >> ')

    # time.sleep(3)
    # plt.close()

    if len(fibers) > 1:
        print('FEU Position for maximum transmission in all fibers:', np.round(xpts[np.argmax(transmissions)],3))
        return np.round(xpts[np.argmax(transmissions)],4), np.nanmax(transmissions)
    else:
        return np.round(xpts[np.argmax(int_flux)],4)

def slit_scan(start=4.32, stop=4.40, step=0.005, plot=True, fibers=None, plotTrace=False,axis='y', N_order_local=7):
    '''
    Scan over given FEU TTM region to maximize flux
    '''
    Path = Organization.get_path('Spec')
    Path = list(Path)
    Path[1] =Path[1].split('.')[0]
    ### Take backgrounds
    bkgd, smoothed_thermal_noise, badpixmap = _take_backgrounds(itime, coadds, nframes=5)
    badpixmap = np.ones((2048,2048))
    # badpixcube = np.expand_dims(badpixmap, axis=0)
    # bkgd = np.zeros((2048, 2048))
    # smoothed_thermal_noise = np.zeros((2048, 2048))

    ### Set up scanning points
    xpts, datacube = _take_slit_scan(start=start, stop=stop, step=step, bkgd=bkgd, Path=Path, fibers=fibers,axis=axis)

    best_pos = _analyze_slit_scan(datacube=datacube, xpts=xpts, badpixmap=badpixmap, plot=plot, Path=Path, fibers=fibers, plotTrace=plotTrace, N_order_local=N_order_local)
    
    ## take new frames
    #new_datacube = np.empty((3, 2048, 2048))
    #spec.TakeImage(Coadds=coadds, tint=itime)
    #new_datacube[i] = spec.im[0].T[:,::-1] - bkgd
    
    return best_pos

def _take_fam_scan(start=-3, stop=3, step=1, fibers=None, bkgd=None, Path=None, nreads=10, verbose=True):
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
        Path = ('/nfiudata/FAM_Scans/', '')
    Track_Flag = track.is_tracking()
    ### This isn't working
    if Track_Flag:
        goal0 = int(track.get_goal()[0])
    try:
        if bkgd is None:
            bkgd = _take_pd_backgrounds(nreads=nreads)
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
                if not track.is_tracking(): 
                    track.start_tracking()
                _ = _acquire_fiber(fibr, verbose=verbose)
                print('tracking again...')
        
                # track.stop_tracking()
            # track.start_tracking()
            time.sleep(1)
            ### Use average position for scan center
            ipos = np.empty((2,10))
            for i in range(10):
                pos = FAM.get_pos()
                #print(i, pos)
                ipos[0,i] = pos[0]
                ipos[1,i] = pos[1]
                time.sleep(0.1)
            ipos = np.mean(ipos, axis=-1)
            # print('mean FAM pos:')
            # print(ipos)

            xpixtorad = FAM.pix2ttmUx
            ypixtorad = FAM.pix2ttmUy
            xpts = xpixtorad*np.arange(start, stop+step/2., step)+ipos[0]
            ypts = ypixtorad*np.arange(start, stop+step/2., step)+ipos[1]
            
            track.stop_tracking()
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
                        time.sleep(0.01)
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
            
            # print('')
            # print('')
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
        Path = ('/nfiudata/FAM_Scans/', '')
    if fibers is None:
        fibers = [track.get_goal()[0]]
    
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
        # plt.xlabel('PSF offset X direction (pixel)')
        # plt.ylabel('PSF offset Y direction (pixel)')
        plt.xlabel('PSF offset Y direction (pixel)')
        plt.ylabel('PSF offset X direction (pixel)')
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
        # plt.xlabel('PSF offset X direction (pixel)')
        # plt.ylabel('PSF offset Y direction (pixel)')
        plt.xlabel('PSF offset Y direction (pixel)')
        plt.ylabel('PSF offset X direction (pixel)')
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
        #plt.xlabel('PSF offset X direction (pixel)')
        #plt.ylabel('PSF offset Y direction (pixel)')
        plt.xlabel('PSF offset Y direction (pixel)')
        plt.ylabel('PSF offset X direction (pixel)')
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
        # plt.xlabel('PSF offset X direction (pixel)')
        # plt.ylabel('PSF offset Y direction (pixel)')
        plt.xlabel('PSF offset Y direction (pixel)')
        plt.ylabel('PSF offset X direction (pixel)')
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
        psfpars = track.get_psf_cent()
        xvals[i] = psfpars[0]
        yvals[i] = psfpars[1]
        time.sleep(0.02)
    xpix = np.round(np.mean(xvals),2)
    ypix = np.round(np.mean(yvals),2)

    track.set_fiber_loc((xpix, ypix), int(fnum))
    track._loadGoals()
    print('Updated fiber position:', track.goals[int(fnum)])

def fam_scan(start=-4, stop=4, step=1, fibers=None, savename='', bkgd=None, plot=True, fit=True, overwrite=False, update_fiber=False, nreads=10):
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
    Path = ('/nfiudata/FAM_Scans/',tstamp)
    ### Take backgrounds
    # modified March 11, 2023. Put into fiber_finding()
    # bkgd = _take_pd_backgrounds(nreads=nreads)
    
    xpts, ypts, scans, cred2fluxes = _take_fam_scan(start=start, stop=stop, step=step, fibers=fibers, bkgd=bkgd, Path=Path, nreads=nreads)

    xoffset, yoffset = _analyze_fam_scan(scan=scans, xpts=xpts, ypts=ypts, cred2fluxes=cred2fluxes, 
                                plot=plot, fit=fit, Path=Path, savename=savename, fibers=fibers)

    if update_fiber:
        newx = np.round(xoffset[0],3)
        newy = np.round(yoffset[0],3)
        fnum = track.get_psf_goal()[-1]
        ### Make sure tracking is off while we reset fiber positions
        Track_Flag = track.is_tracking()
        if Track_Flag:
            track.stop_tracking()
        confirm = input('Update position of fiber '+str(int(fnum))+' to '+str(newx)+', '+str(newy)+'? [Y/n] >> ')
        if confirm in ['Y', 'Yes', 'yes']:
            _update_fiber_position(fnum, (xoffset[0], yoffset[0]))
        else:
            print('Fiber positions not updated')
            print((xoffset[0], yoffset[0]))
        if Track_Flag:
            track.start_tracking()

    return xoffset, yoffset

def scan_mode(noll, start=-0.03, stop=0.03, step=0.01, grid_size=13, stepsize=2, refsrf=nominalFlat, fiber=2, nreads=10, bkgd=None,
             show=True, date_dir='', verbose=True):
    ''' Run scan over sepcified zernike mode
    
    Arguments:
    noll (int) Zernike noll to scan over

    Keyword arguments:
    start (float, default -0.03): Lower limit of aberattion axis. Units are whatever DM.py uses, IDK
    stop (float, default 0.03): Upper limit on aberration axis (inclusive)
    step (float, default 0.01): Stepsize along aberration axis
    gridsize (int, default 13): Size of TT scan do perform, in steps. Default gives a 13x13
    stepsize (float, default 2): Stepsize of the TT scan, in cred 2 pixels
    refsrf (string, default nominalFlat) Name of reference DM surface to apply aberration to. If
            'current', uses whatever the DM shape is at the start of the scan (useful for iterative tuning)

    Outputs: Tuple of scan x points, y points, applied amplitudes, pd readings, and cred2 fluxes. Path is
    /home/nfiudev/dev/FAM_Scans/Noll_<noll>_PD_Values.npz. Note that this will overwrite any preivous scan of that noll!!!

    Returns: Tuple of scan x points, y points, applied amplitudes, pd readings, and cred2 fluxes.
    '''
    ### Reset to initial DM shape
    # K2AO.load_dm_shape(dmfname)
    # inimap = K2AO.DM_shape_map
    ### Loop through points
    startTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())               
    tstamp = startTime.replace('-','')                                          
    tstamp = tstamp.replace(':','')                                             
    tstamp = tstamp.replace(' ','')                                             
    tstamp = tstamp[2:-2]
    
    try:
        ### Get DM reference. Either a flatmap or 'current', do do everything relative to the current
        ### DM setting
        if refsrf != 'current':
            flat = DM.setFlatSurf(refsrf)
        else:
            flat = DM.getSurf()

        # Use bkgd from before the scan. Loaded as argument.
        # bkgd = _take_pd_backgrounds(nreads=nreads)
        # time.sleep(0.5)
        zpts = np.arange(start, stop+step/2., step)
        scans = np.empty((zpts.size, grid_size, grid_size))
        cred2s = np.empty((zpts.size, grid_size, grid_size))
        for i in range(zpts.size):
            if verbose:
                print('Amplitude:', np.round(zpts[i],3))
            ### Apply zernike
            dmmap = DM.pokeZernike(zpts[i], noll, bias=flat)
            ### Read photometer
            xpts, ypts, scans[i], cred2s[i] = _take_fam_scan(bkgd=bkgd, start=-int((stepsize*(grid_size-1))/2), stop=int((stepsize*(grid_size-1))/2), step=stepsize, nreads=nreads, verbose=verbose)
            ### Set DM back where it started
            DM.setSurf(flat)
            _acquire_fiber(fiber, toggle_track=True)

        ### Save data
        Path = os.path.join('/home/nfiudev/dev/FAM_Scans', date_dir)

        fullname = os.path.join(Path, 'Noll_%02d_PD_Values.npz' %(noll))
        print('Saving to:', fullname)
        np.savez(fullname, xpts=xpts, ypts=ypts, amps=zpts, pdfluxes=scans, crdfluxes=cred2s)
        
        # make the plot
        amp = plot_zernike_scan_peaks(noll, show=show, date_dir=date_dir)
        rounded_amp = np.round(amp, 6)
        print('Noll ' + str(noll) + ': amp = ' +str(rounded_amp))
        
        ### return
        # return xpts, ypts, zpts, scans, cred2s
        return rounded_amp

    ### Reload inital DM shape if something went wrong
    except(Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        DM.setSurf(flat)
        print('DM shape reset')
        raise e

def plot_zernike_scan_peaks(noll, fit=True, show=True, date_dir = ''):
    ''' Plots raw peak voltage of scan vs applied amplitude

    Arguments: noll (int): noll index to load. Assume default filename used by scan_mode()

    Keyword Arguments: fit (bool, default True): Whether to fit a quadratic
    '''
    
    # Path = os.path.join('/home/nfiudev/dev/FAM_Scans', date_dir)
    Path = os.path.join('/nfiudata/FAM_Scans', date_dir)
    
    startTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
    tstamp = startTime.replace('-','')
    tstamp = tstamp.replace(':','')
    tstamp = tstamp.replace(' ','')
    tstamp = tstamp[2:-2]
    fullname = os.path.join(Path, 'Noll_%02d_PD_Values.npz' %(noll))
    scandata = np.load(fullname)
    xpts, ypts = scandata['xpts'], scandata['ypts']
    amps = scandata['amps']
    pdfluxes = scandata['pdfluxes']
    crdfluxes = scandata['crdfluxes']

    peaks = np.empty(pdfluxes.shape[0])

    for i in range(peaks.size):
        peaks[i] = np.max(pdfluxes[i])

    fig, ax = plt.subplots()
    ax.scatter(amps, peaks, color='k')
    if fit:
        mpts = np.linspace(np.nanmin(amps), np.nanmax(amps), 1000)
        gauss = models.Gaussian1D(amplitude=(np.nanmax(peaks)-np.nanmin(peaks)), mean=0, stddev=0.1)
        fitter = modeling.fitting.LevMarLSQFitter()
        fitted = fitter(gauss, amps[np.isfinite(peaks)], peaks[np.isfinite(peaks)])
        fitcurve = fitted(mpts)
        ax.plot(mpts, fitcurve, color='k', linestyle='--')
        bestpt = mpts[np.argmax(fitcurve)]
        ax.axvline(bestpt, color='r', linestyle='--')
        ax.text(bestpt+0.015, np.min(peaks)+0.5*(np.max(peaks)-np.min(peaks)), 'Fit peak = '+str(np.round(bestpt,3))+' V')
    np.savez( os.path.join(Path,tstamp+'_noll'+str(noll)+'_rawpeaks.npz'), amps=amps, peaks=peaks)
    ax.set_title(nolls[noll] + ' Raw Voltage')
    ax.set_xlabel('Mode Amplitude')
    ax.set_ylabel('Peak PD voltage')
    plt.savefig(os.path.join(Path,tstamp+'_noll'+str(noll)+'_rawpeaks.png'), bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close('all')

    if fit: return bestpt
    return amps[np.argmax(peaks)]

def plot_zernike_scans(noll):
    '''Plots 2D images of each frame in scan

    Arguments: noll (int): noll index to load. Assume default filename used by scan_mode()
    '''
    #Path = '/home/nfiudev/dev/FAM_Scans/'
    Path = '/nfiudata/FAM_Scans/'
    startTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
    tstamp = startTime.replace('-','')
    tstamp = tstamp.replace(':','')
    tstamp = tstamp.replace(' ','')
    tstamp = tstamp[2:-2]
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
    plt.savefig(Path+tstamp+'_Noll_'+str(noll).zfill(2)+'_raw_TT_scans.png',bbox_inches='tight')
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
    Path = '/nfiudata/FAM_Scans/'
    startTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
    tstamp = startTime.replace('-','')
    tstamp = tstamp.replace(':','')
    tstamp = tstamp.replace(' ','')
    tstamp = tstamp[2:-2]
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
        mpts = np.linspace(np.nanmin(amps), np.nanmax(amps), 1000)
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
    np.savez(Path+tstamp+'_noll'+str(noll)+'_fitpeaks.npz', amps=amps, peaks=peaks)
    plt.title(nolls[noll]+' Best-fit Amplitude')
    plt.xlabel('Mode Amplitude')
    plt.ylabel('Best-Fit Peak PD voltage')
    plt.savefig(Path+tstamp+'_noll'+str(noll)+'_fitpeaks.png', bbox_inches='tight')
    plt.show()

    if fit: return bestpt
    return amps[np.argmax(peaks)]

def set_dm_shape(noll, zpt, refsrf=nominalFlat):
    ''' Set DM shape to specified noll/amplitude

    Arguments:
    noll (int): noll index of zernike to apply
    zpt (float): amplitude of aberration to apply, in DM.py units

    Keyword arguments:
    refsrf (str, default nominalFlat): DM reference surface to apply
            aberration on top of. If 'current', uses the current shape
            of the DM (useful for iterative tuning)
    '''
    ### Reset to initial DM shape
    if refsrf != 'current':
        flat = DM.setFlatSurf(nominalFlat)
    else:
        flat = np.copy(DM.getSurf())
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
        DM.setFlatSurf(nominalFlat)
        print('DM shape reset')
        raise e


def reset_dm_shape():
    '''Resets the DM shape to the 220225_601 flat map, in case you fuck up
    '''
    DM.setFlatSurf(nominalFlat)

def set_kpic_laser(state):
    if state == 'on':
        print('Turning on 2 micron laser...')
        tools.laser_connect()
        tools.laser_setcurrent(250)
        tools.laser_enable()
    elif state == 'off':
        print('Turning off 2 micron laser...')
        tools.laser_connect()
        tools.laser_setcurrent(0)
        tools.laser_disable()

### Scripts to run calibs as pipeline ###

# this needs a lot of checking and work. Do SAM rot and PSM focus have connect commands??
def check_kpic_devices():
    # connect all 12 devices
    if not SAM.is_connected():
        SAM.connect()
    
    if not PSM.is_connected():
        PSM.connect()
        while not PSM.is_connected():
            time.sleep(0.5)
        PSM.close_loop()

    if not FAM.is_connected():
        FAM.connect()
        while not FAM.is_connected():
            time.sleep(0.5)
        FAM.close_loop()

    if not modeChange.is_connected():
        modeChange.connect()
    if not LightSrcRet.is_connected():
        LightSrcRet.connect()
    if not cam.is_connected():
        cam.connect()
    if not Coronagraph.is_connected():
        Coronagraph.connect()
        while not Coronagraph.is_connected():
            time.sleep(0.5)
        Coronagraph.close_loop()
    if not PIAA.is_connected():
        PIAA.connect()
        # PIAA is in open loop by default
    if not PyWFSPickoff.is_connected():
        PyWFSPickoff.connect()
    # if not PyWFSFilt.is_connected():
    #     PyWFSFilt.connect()
    #     while not PyWFSFilt.is_connected():
    #         time.sleep(0.5)
    #     PyWFSFilt.close_loop()
    if not FilterWh.is_connected():
        FilterWh.connect()
    if not FIU_Fiber.is_connected():
        FIU_Fiber.connect()
        while not FIU_Fiber.is_connected():
            time.sleep(0.5)
        FIU_Fiber.close_loop()

    # check control script is active
    if (SAM.is_rot_active() and SAM.is_active()):
        pass
    else:
        SAM.activate_control_script()
    if (PSM.is_foc_active() and PSM.is_active()):
        pass
    else:
        PSM.activate_control_script()
    if FAM.is_active():
         pass
    else:
        FAM.activate_control_script()
    if LightSrcRet.is_active():
        pass
    else:
        LightSrcRet.activate_control_script()
    if modeChange.is_active():
         pass
    else:
        modeChange.activate_control_script()
    if cam.is_active():
         pass
    else:
        cam.activate_control_script()
    if Coronagraph.is_active():
         pass
    else:
        Coronagraph.activate_control_script()
    if PIAA.is_active():
         pass
    else:
        PIAA.activate_control_script()
    if PyWFSFilt.is_active():
         pass
    else:
        PyWFSFilt.activate_control_script()
    if PyWFSPickoff.is_active():
         pass
    else:
        PyWFSPickoff.activate_control_script()
    if FilterWh.is_active():
         pass
    else:
        FilterWh.activate_control_script()
    if FIU_Fiber.is_active():
         pass
    else:
        FIU_Fiber.activate_control_script()

    time.sleep(10)
    print('All devices connected and active.')

# Run setup for calibrations. Starts with internal KPIC source
def calib_setup(pywfsd, light_source='kpic', psm_foc=2.54, k2ao_pos='open'):
    ''' 
    Initial setup, step 1 of the calibration document
    By default turns on the kpic light sources (laser + MIR lamp)
    pywfsd: PyWFS dichroic position (ds or empty1)
    '''

    # First check that all stages we need are active or in closed loop
    check_kpic_devices()

    # SAM rotator in, and SAM to pd sf1 first
    if not SAM.get_rot_named_pos() == 'in':
        SAM.set_rot_pos('in', block=True)

    # make sure SAM rot is in
    # if SAM.get_rot_named_pos() == 'in':
    #     SAM.set_pos('pd_sf1', block=True)

    # PSM focus and bundle
    if (np.round(PSM.get_foc_pos(), 3)[0] - psm_foc) > 1e-2:
        print('Setting PSM focus to ' + str(psm_foc))
        PSM.set_foc_pos(psm_foc, foc=2, block=True)
    
    # Both filter wheels to h
    if FilterWh.get_named_pos() != 'h':
        FilterWh.set_pos('h')
    if PyWFSFilt.get_named_pos() != 'h':
        PyWFSFilt.set_pos('h')

    # Check we are in focal viewing mode
    if not modeChange.get_named_pos() == 'focus':
        modeChange.set_pos('focus', block=True)

    # DS or empty1 (no dichroic)
    if pywfsd == 'ds' or pywfsd == 'empty1':
        if not PyWFSPickoff.get_named_pos() == pywfsd:
            PyWFSPickoff.set_pos(pywfsd, block=True)
    else:
        print('Check input pyWFSd. It should be ds or empty1.')

    if PyWFSPickoff.get_named_pos() == 'ds':
        if not Coronagraph.get_named_pos() == 'pupil_mask':
            Coronagraph.set_pos('pupil_mask')

    elif PyWFSPickoff.get_named_pos() == 'empty1':
        if not Coronagraph.get_named_pos() == 'pypo_out':
            Coronagraph.set_pos('pypo_out')

        # corona_pos = Coronagraph.get_pos()
        # # if far from goal position
        # if (corona_pos[0] - 3.3 > 1e-3 ) or (corona_pos[1] - 1.88 > 1e-3 ):
        #     Coronagraph.set_pos((3.3, 1.88))
        #     print('Set coronagraph to (3.3, 1.88).')

    # Turn on PD Pickoff
    power_on_pdp = 'modify -s k2aopower OUTLET_HC7=on'
    os.system(power_on_pdp)
    time.sleep(2)
    flip_to_pd = 'modify -s kpic PDPONAM=PD'
    os.system(flip_to_pd)
    time.sleep(1)

    if light_source == 'kpic':
        # light Src retractor in
        if not LightSrcRet.get_named_pos() == 'in':
            LightSrcRet.set_pos('in', block=False)
        # Turn on laser
        set_kpic_laser('on')
        # Turn on MIR Lamp (usually done with pdu_gui)
        power_on_mir = 'modify -s k2aopower OUTLET_HB3=on'
        os.system(power_on_mir)

        # put on flat map for kpic source
        internal_flat = DM.setFlatSurf(nominalFlat)
        print('Applied DM map ' + nominalFlat + ' for internal source.')

    # Need to be careful to make sure we have K2AO bench access
    elif light_source == 'k2ao':
        # light Src retractor in
        if not LightSrcRet.get_named_pos() == 'out':
            LightSrcRet.set_pos('out', block=False)

        k2ao_access = input('Confirm access to K2AO bench to turn on cal source? (Y/n) >>> ')
        if k2ao_access == 'Y':
            print('Turning on K2AO source is on, with position=' + k2ao_pos)
            ktl.write('ao2', 'obswon', '1')
            Light_Src_cmds.set(k2ao_pos)
        else:
            print('Not turning on K2AO SFP.')

    # start with PSF in center position
    FAM.set_pos('center', block=True)

    # zero any offsets that might be there
    # track.set_astrometry(0, 0)
    track.reset_astrometry()
    time.sleep(1)
    track.use_disp(False)  # no DAR

def calib_setup_sfp(pywfsd, light_source='kpic_sfp', psm_foc=2.54, k2ao_pos='open', 
move_fsm=False, fsm_x=6.3037, fsm_y=6.8104):
    ''' 
    Initial setup, step 1 of the calibration document
    By default turns on the kpic light sources (laser + MIR lamp) and puts them
    in front of K2AO.

    pywfsd: PyWFS dichroic position (ds or empty1)
    '''

    # First check that all stages we need are active or in closed loop
    check_kpic_devices()

    # SAM rotator in, and SAM to pd sf1 first
    if not SAM.get_rot_named_pos() == 'in':
        SAM.set_rot_pos('in', block=True)

    # make sure SAM rot is in
    # if SAM.get_rot_named_pos() == 'in':
    #     SAM.set_pos('pd_sf1', block=True)

    # PSM focus and bundle
    if (np.round(PSM.get_foc_pos(), 3)[0] - psm_foc) > 1e-2:
        print('Setting PSM focus to ' + str(psm_foc))
        PSM.set_foc_pos(psm_foc, foc=2, block=True)
    
    # Both filter wheels to h
    if FilterWh.get_named_pos() != 'h':
        FilterWh.set_pos('h')
    if PyWFSFilt.get_named_pos() != 'h':
        PyWFSFilt.set_pos('h')

    # Check we are in focal viewing mode
    if not modeChange.get_named_pos() == 'focus':
        modeChange.set_pos('focus', block=True)

    # DS or empty1 (no dichroic)
    if pywfsd == 'ds' or pywfsd == 'empty1':
        if not PyWFSPickoff.get_named_pos() == pywfsd:
            PyWFSPickoff.set_pos(pywfsd, block=True)
    else:
        print('Check input pyWFSd. It should be ds or empty1.')

    if PyWFSPickoff.get_named_pos() == 'ds':
        if not Coronagraph.get_named_pos() == 'pupil_mask':
            Coronagraph.set_pos('pupil_mask')

    elif PyWFSPickoff.get_named_pos() == 'empty1':
        if not Coronagraph.get_named_pos() == 'pypo_out':
            Coronagraph.set_pos('pypo_out')

        # corona_pos = Coronagraph.get_pos()
        # # if far from goal position
        # if (corona_pos[0] - 3.3 > 1e-3 ) or (corona_pos[1] - 1.88 > 1e-3 ):
        #     Coronagraph.set_pos((3.3, 1.88))
        #     print('Set coronagraph to (3.3, 1.88).')

    # Turn on PD Pickoff
    power_on_pdp = 'modify -s k2aopower OUTLET_HC7=on'
    os.system(power_on_pdp)
    time.sleep(2)
    flip_to_pd = 'modify -s kpic PDPONAM=PD'
    os.system(flip_to_pd)
    time.sleep(1)

    if light_source == 'kpic_int':
        # light Src retractor in
        if not LightSrcRet.get_named_pos() == 'in':
            LightSrcRet.set_pos('in', block=False)
        # Turn on laser
        set_kpic_laser('on')
        # Turn on MIR Lamp (usually done with pdu_gui)
        power_on_mir = 'modify -s k2aopower OUTLET_HB3=on'
        os.system(power_on_mir)
        
        show_swit = 'show -s kpic LSRCSWIT'
        print('Current switch position: ')
        swit_curpos = sp.getoutput(show_swit)
        print(swit_curpos)
        if not 'KPICIN' in swit_curpos:
            print('Not in KPICIN position. Moving to KPICIN...')
            move_to_kpicin = 'modify -s kpic LSRCSWIT=KPICIN'
            os.system(move_to_kpicin)
    
        # put on flat map for kpic source
        internal_flat = DM.setFlatSurf(nominalFlat)
        print('Applied DM map ' + nominalFlat + ' for internal source.')

    # Need to be careful to make sure we have K2AO bench access
    elif light_source == 'kpic_sfp':
        # light Src retractor in
       # if not LightSrcRet.get_named_pos() == 'out':
        #    LightSrcRet.set_pos('out', block=False)

        k2ao_access = input('Confirm access to K2AO bench to put KPIC sources in front of K2AO? (Y/n) >>> ')
        if k2ao_access == 'Y':
            print('Turning on KPIC sources and putting them in front of K2AO.')
           # ktl.write('ao2', 'obswon', '1')
            #Light_Src_cmds.set(k2ao_pos)
            kpic_sfp_source(source='on', move_fsm=move_fsm, fsm_x=fsm_x, fsm_y=fsm_y)
        else:
            print('Not moving KPIC sources in front of K2AO.')

    # start with PSF in center position
    #FAM.set_pos('center', block=True)

    # zero any offsets that might be there
    # track.set_astrometry(0, 0)
    # track.reset_astrometry()
    # time.sleep(1)
    # track.use_disp(False)  # no DAR

def kpic_offsky(cal_sources_off=True, sam_rot_out=False, zero_dm=True):
    ''' 
    Default off sky positions. To run after calibrations and before going on-sky
    pywfsd: ds or empty1
    '''

    # First check that all stages we need are active or in closed loop
    # check_kpic_devices()
    
    # turn off tracking and DAR, and remove offset
    track.stop_tracking()
    time.sleep(1)
    # track.set_astrometry(0, 0)
    track.reset_astrometry()
    time.sleep(1)
    track.use_disp(False)

    # light Src retractor in
    if LightSrcRet.get_named_pos() == 'in':
        LightSrcRet.set_pos('out', block=True)

    # SAM rotator out. This is needed for regular NIRSPAO to have light
    # but SAM rot must be in for good NIRSPEC backgrounds (otherwise fibers are not in slit)
    if sam_rot_out:
        if SAM.get_rot_named_pos() == 'in':
            SAM.set_rot_pos('out', block=True)

    # Check we are in focal viewing mode
    if not modeChange.get_named_pos() == 'focus':
        modeChange.set_pos('focus', block=True)
    
    # option to run without dichroic (only possible for SH WFS mode)
    # ds for dichroic in, empty1 for dichroic out
    # if not PyWFSPickoff.get_named_pos() == pywfsd:
    #     PyWFSPickoff.set_pos(pywfsd, block=True)

    # Both filter wheels to h
    if FilterWh.get_named_pos() != 'h':
        FilterWh.set_pos('h')
    if PyWFSFilt.get_named_pos() != 'h':
        PyWFSFilt.set_pos('h')

    # turn off PD Pickoff - should be off already
    # flip_to_nirspec = 'modify -s kpic PDPONAM=NIRSPEC'
    # os.system(flip_to_nirspec)
    # time.sleep(2)
    if cal_sources_off:
        power_off_pdp = 'modify -s k2aopower OUTLET_HC7=off'
        os.system(power_off_pdp)

        # Turn off kpic laser
        set_kpic_laser('off')
        # Turn off MIR Lamp (usually done with pdu_gui)
        power_off_mir = 'modify -s k2aopower OUTLET_HB3=off'
        os.system(power_off_mir)

        # turn off k2ao source
        # only if we have access to full K2AO bench, so an option
        # if keck_source_off:
        # k2ao_access = input('Confirm access to K2AO bench to turn off cal source? (Y/n) >>> ')
        # if k2ao_access == 'Y':
        ktl.write('ao2', 'obswon', '0')
        # else:
        #     print('Not turning off K2AO SFP.')

    # Zero the DM
    if zero_dm:
        DM.zeroAll()
        print('DM voltages zeroed.')
    print('KPIC is in off-sky state. Remember to take NIRSPEC backgrounds if you have time.')

def kpic_onsky(dmmap_file):
    ''' 
    Default on sky positions. To run after calibrations and before going on-sky
    dmmap_file: map to put on DM. Determined from daycals.
    '''

    # turn off tracking, remove offset to begin with. But turn on DAR.
    track.stop_tracking()
    time.sleep(1)
    # track.set_astrometry(0, 0)
    track.reset_astrometry()
    time.sleep(1)

    dar_ready = input('Ready to turn on DAR? (Y/n)')
    if dar_ready == 'Y':
        track.use_disp(True)

    # light Src retractor out
    if LightSrcRet.get_named_pos() == 'in':
        LightSrcRet.set_pos('out', block=True)

    # SAM rotator back in.
    if SAM.get_rot_named_pos() == 'out':
        SAM.set_rot_pos('in', block=True)

    # Check we are in focal viewing mode
    if not modeChange.get_named_pos() == 'focus':
        modeChange.set_pos('focus', block=False)

    # probably don't touch this
    # if not PyWFSPickoff.get_named_pos() == pywfsd:
    #     PyWFSPickoff.set_pos(pywfsd)
    
    # Filter wheel to h
    if FilterWh.get_named_pos() != 'h':
        FilterWh.set_pos('h', block=False)

    # Turn off PD Pickoff (should be in NIRSPEC position)
    power_off_pdp = 'modify -s k2aopower OUTLET_HC7=off'
    os.system(power_off_pdp)

    # Turn off kpic laser
    set_kpic_laser('off')
    # Turn off MIR Lamp (usually done with pdu_gui)
    power_off_mir = 'modify -s k2aopower OUTLET_HB3=off'
    os.system(power_off_mir)
    # turn off k2ao source
    ktl.write('ao2', 'obswon', '0')
    
    # Set the DM map
 #   flat_ds = np.load(dmmap_file)
    flat_ds = np.load(os.path.join(DM.flatdir, dmmap_file)) 
    DM.setSurf(flat_ds)
    print('KPIC is ready for on-sky observations. Applied DM map ' + dmmap_file)

def use_keck_source(k2ao_pos='open', move_fsm=False, fsm_x = 6.32451, fsm_y = 6.83247):

    # first turn off tracking
    track.stop_tracking()
    time.sleep(1)
    # track.set_astrometry(0, 0)
    track.reset_astrometry()

    # move DFB to mirror
    show_dfb = 'show -s ao obdbname'
    print('Current DFB position: ')
    dfb_curpos = sp.getoutput(show_dfb)
    print(dfb_curpos)
    if not 'mirror' in dfb_curpos:
        print('Not in mirror position. Moving to mirror...')
        move_to_mirror = 'modify -s ao obdbname=mirror'
        os.system(move_to_mirror)

    # TO-DO: check if kpic sources are on
    # Turn off kpic laser
    set_kpic_laser('off')
    # Turn off MIR Lamp (usually done with pdu_gui)
    power_off_mir = 'modify -s k2aopower OUTLET_HB3=off'
    os.system(power_off_mir)

    # light source retractor out
    if LightSrcRet.get_named_pos() == 'in':
        LightSrcRet.set_pos('out', block=True)

    # turn on k2ao source
    ktl.write('ao2', 'obswon', '1')
    Light_Src_cmds.set(k2ao_pos)
    print('K2AO source is on, with ' + k2ao_pos)

    # start with PSF in center position - important to align FSM2 below
    FAM.set_pos('center', block=True)

    # move FSM2 (pyWFS FSM) to desired position
    if move_fsm:
        execute_fsm_move(fsm_x, fsm_y)


def kpic_sfp_source(source='on', move_fsm=False, fsm_x = 6.3037, fsm_y = 6.8104):

    # first close NIRC2 shutters
    show_shutter = nirc2.get_shutter()
    print(f'Current shutter position: {show_shutter}')
    if not 'closed' in nirc2.get_shutter():
        print('NIRC2 shutter not in closed position. Closing shutter...')
        nirc2.close_shutter()

    # then turn off tracking
    track.stop_tracking()
    time.sleep(1)
    # track.set_astrometry(0, 0)
    track.reset_astrometry()

    # move DFB to mirror
    show_dfb = 'show -s ao obdbname'
    print('Current DFB position: ')
    dfb_curpos = sp.getoutput(show_dfb)
    print(dfb_curpos)
    if not 'mirror' in dfb_curpos:
        print('Not in mirror position. Moving to mirror...')
        move_to_mirror = 'modify -s ao obdbname=mirror'
        os.system(move_to_mirror)

    # Toggle optical switch
    show_swit = 'show -s kpic LSRCSWIT'
    print('Current switch position: ')
    swit_curpos = sp.getoutput(show_swit)
    print(swit_curpos)
    if source == 'off':
        if not 'KPICIN' in swit_curpos:
            print('Not in KPICIN position. Moving to KPICIN...')
            move_to_kpicin = 'modify -s kpic LSRCSWIT=KPICIN'
            os.system(move_to_kpicin)
    else:
        if not 'K2AOIN' in swit_curpos:
            print('Not in K2AOIN position. Moving to K2AOIN...')
            move_to_k2aoin = 'modify -s kpic LSRCSWIT=K2AOIN'
            os.system(move_to_k2aoin)

    if source == 'off':
        print('Turning off laser and MIR lamp...')
        # Turn off kpic laser
        set_kpic_laser('off')
        # Turn off MIR Lamp (usually done with pdu_gui)
        power_off_mir = 'modify -s k2aopower OUTLET_HB3=off'
        os.system(power_off_mir)

    # light source retractor out
    if LightSrcRet.get_named_pos() == 'in':
        print('Moving light source retractor out...')
        LightSrcRet.set_pos('out', block=True)

    if source == 'on': 
        # start with PSF in center position - important to align FSM2 below
        print('Setting FAM position to center...')
        FAM.set_pos('center', block=True)

        print('Turning on laser and MIR lamp...')
         # Turn on kpic laser
        set_kpic_laser('on')
        # Turn off MIR Lamp (usually done with pdu_gui)
        power_on_mir = 'modify -s k2aopower OUTLET_HB3=on'
        os.system(power_on_mir)


    if source == 'off':
        # turn on k2ao source
        ktl.write('ao2', 'obswon', '1')
        print('K2AO source is on')
    else:
        # turn off k2ao source
        k2ao_pos = 'open'
       # ktl.write('ao2', 'obswon', '1')
        Light_Src_cmds.set(k2ao_pos)
        ktl.write('ao2', 'obswon', '0')
        print('K2AO source is off, with ' + k2ao_pos)

    if source == 'off':
        # Move SFP pos to default pws pos
        print('Moving SFP position to pws...')
        sfp.set_name_position('pws')
    else:
        # Move SFP position to kpiccal
        print('Moving SFP position to kpiccal...')
        sfp.set_name_position('kpiccal')

    # move FSM2 (pyWFS FSM) to desired position
    if move_fsm:
        execute_fsm_move(fsm_x, fsm_y)



# Doesn't work perfectly. Need to debug
def execute_fsm_move(fsm_x, fsm_y):
    '''
    Move FSM2 (PyWFS FSM) to desired x and y positions
    '''
    _x_pos = sp.getoutput('show -s ao pmcf2x')
    cur_x_pos = float(_x_pos.strip().split(' mm')[0].split('= ')[1])
    _y_pos = sp.getoutput('show -s ao pmcf2y')
    cur_y_pos = float(_y_pos.strip().split(' mm')[0].split('= ')[1])

    print('Current x, current y:')
    print(cur_x_pos, cur_y_pos)

    # compute the deltas to move by
    _delta_x = fsm_x - cur_x_pos
    _delta_y = fsm_y - cur_y_pos

    # convert to strings and apply one at a time
    delta_x = str(_delta_x)
    delta_y = str(_delta_y)

    print(_delta_x, _delta_y)

    # only apply if meaningful
    # need to take absolute value!! March 11, 2023
    if np.abs(_delta_x) > 0.0005:
        # x move
        os.system('modify -s ao pmcf2rx='+delta_x)
        time.sleep(0.5)
        os.system('modify -s ao pmcf2ry=0')
        # time.sleep(0.5)
        fsm_ready = input('Apply FSM move of '+delta_x+' in x? (Y/n) >>> ')
        if fsm_ready == 'Y':
            os.system('modify -s ao pmcf2gr=1')
        else:
            print('Not moving in x.')
        time.sleep(2)

    if np.abs(_delta_y) > 0.0005:
        # y move
        os.system('modify -s ao pmcf2ry='+delta_y)
        time.sleep(0.5)
        os.system('modify -s ao pmcf2rx=0')
        fsm_ready = input('Apply FSM move '+delta_y+' in y? (Y/n) >>> ')
        if fsm_ready == 'Y':
            os.system('modify -s ao pmcf2gr=1')
        else:
            print('Not moving in y.')
        time.sleep(2)

    final_x_pos = sp.getoutput('show -s ao pmcf2x')
    final_y_pos = sp.getoutput('show -s ao pmcf2y')

    print('Current FSM positions (x, y):')
    print(final_x_pos.strip(), final_y_pos.strip())

def setup_for_scam(ini_foc=2.52, scam_nd = 'nd6', Scam_tint=1500, scan_fib=2):

    # move to nd6 (default) to prepare for PSM scan on scam
    Light_Src_cmds.set(scam_nd)
    print('Inserted ' + scam_nd)
    time.sleep(1)

    # increase cred2 tint now to keep PSF visible
    # CHECK THIS
    cam.set_fps(2)
    time.sleep(1)
    cam.set_tint(0.497)  # this is equal to Hmag=9.2 in GUI
    print('Increased CRED2 exp time, equivalent of Hmag=9.2 preset.')
    print('Taking CRED2 background')
    
    # take a CRED2 background
    track.stop_tracking()
    FAM.set_pos(fam_bkgd_pos, block = True)
    time.sleep(1)
    proc.save_dark()
    time.sleep(8)

    FAM.set_pos('center', block=True)
    time.sleep(1)

    # start tracking on SF2
    # if not track.is_tracking():
    track.start_tracking()
    time.sleep(0.5)
    _acquire_fiber(scan_fib)
    time.sleep(3)

    # PD pickoff to NIRSPEC
    flip_to_nirspec = 'modify -s kpic PDPONAM=NIRSPEC'
    os.system(flip_to_nirspec)
    time.sleep(1)

    # load presets 
    SAM.load_presets()
    time.sleep(0.5)
    
    # SAM to preset off_slit position
    SAM.set_pos('off_slit', block=True)
    time.sleep(1)

    # set PSM focus to nominal value
    PSM.set_foc_pos(ini_foc, foc=2, block=True)

    # light source retractor out
    if LightSrcRet.get_named_pos() == 'in':
        LightSrcRet.set_pos('out', block=True)

    # take a scam image
    print('Taking a sample SCAM image. Make sure there is light on SCAM.')
    scam.TakeImage(tint=Scam_tint)

    print('Ready for PSM XY Scan.')

def fiber_finding(fibers=[1,2,3,4], start=-5, stop=5, step=0.5, light_source='k2ao', mode='ds'):
    if light_source == 'kpic':
        internal_flat = DM.setFlatSurf(nominalFlat)
        nreads_fibfind = 10
    else:
        # Start with previous flat map - try to find automatically. Assume they are saved in DM.flatdir
        all_flats = glob.glob(os.path.join(DM.flatdir, 'NCPA_map_' + mode + '*'))
        guess_latest = np.sort(all_flats)[-1]
        print('Suggesting the initial flat map: ' + guess_latest)
        ini_flat_good = input('Is the initial flat map good? (Y/n) >>> ')
        if ini_flat_good == 'Y':
            ini_flat = guess_latest
        else:
            user_ini_flat = input('Please input name of initial flat map >>> ')
            ini_flat = user_ini_flat.strip()
        print('Applying DM Map ' + ini_flat)
        old_flat = np.load(ini_flat)
        DM.setSurf(old_flat)
        nreads_fibfind = 100

    # load presets 
    SAM.load_presets()
    time.sleep(1)

    bkgd = _take_pd_backgrounds(nreads=nreads_fibfind)
    time.sleep(2)

    # Need to start tracking first, otherwise _acquire_fiber hangs
    track.start_tracking()

    for fib in fibers:
        # acquire and track
        _acquire_fiber(fib)
        track.start_tracking()
        # SAM position to pd_sfx
        SAM.set_pos('pd_sf' + str(fib), block=True)
        time.sleep(2)

        fib_pow = tools.read_pd(nreads_fibfind).mean()
        print(fib, np.round(fib_pow, 3))
        if fib_pow < 0.2:
            print(' Power on fiber ' + str(fib) + ' is low. Might want to re-align SAM pd_sf' + str(fib))
            # print('This can be done with tools.sam_2D_scan(start, stop, nsteps)')  # TO-DO: automate this

        # scan
        fam_scan(start=start,stop=stop,step=step,fibers=[fib],update_fiber=True,bkgd=bkgd)

    print('Final goals for the fibers:')
    print(track.goals[1])
    print(track.goals[2])
    print(track.goals[3]) 
    print(track.goals[4])

# Iterate over a range of noll indices
# apply after each scan, and plot results in the end
# also plot PD read over time
def ncpa_scan(start_noll=4, final_noll=10, start=-0.3, stop=0.3, step=0.02, nreads=10, date = '221006', refsrf='current', fiber=2,
                grid_size=1, stepsize=2, verbose=False, mode='ds'):

    # Start with previous flat map - try to find automatically. Assume they are saved in ~/dev/
    all_flats = glob.glob(os.path.join(DM.flatdir, 'NCPA_map_' + mode + '*'))
    guess_latest = np.sort(all_flats)[-1]

    print('Suggesting the initial flat map: ' + guess_latest)
    ini_flat_good = input('Is the initial flat map good? (Y/n) >>> ')
    if ini_flat_good == 'Y':
        ini_flat = guess_latest
    else:
        user_ini_flat = input('Please input name of initial flat map >>> ')
        ini_flat = user_ini_flat.strip()
    print('Applying DM Map ' + ini_flat)

    old_flat = np.load(os.path.join(DM.flatdir, ini_flat))
    DM.setSurf(old_flat)

    # initial pd read of background
    ini_bkgd = _take_pd_backgrounds(nreads=50)
    time.sleep(1)

    Path = '/nfiudata/FAM_Scans/'
    date_dir = os.path.join(Path, date)
    if not os.path.exists(date_dir):
        os.makedirs(date_dir)

    modes_to_scan = np.linspace(start_noll, final_noll, final_noll-start_noll+1, dtype=int)
    all_amps = []
    all_pdread = []

    # start tracking if not already
    if not track.is_tracking():
        track.start_tracking()
        time.sleep(1)
        
    # go to SF2 first
    _acquire_fiber(fiber)
    SAM.set_pos('pd_sf' + str(fiber), block=True)
    time.sleep(2)

    ini_pdread = tools.read_pd(nreads).mean() * 1e4 - ini_bkgd
    all_pdread.append(ini_pdread)
    print('Initial PD read: ' + str(ini_pdread))
    time.sleep(0.5)

    # start scanning modes requested
    for n in modes_to_scan:
        # stop tracking immediately before each scan
        track.stop_tracking()
        this_amp = scan_mode(n, start=start, stop=stop, step=step, grid_size=grid_size, stepsize=stepsize, refsrf=refsrf, fiber=fiber, bkgd=ini_bkgd,
                nreads=nreads, show=True, date_dir=date_dir, verbose=verbose)
        all_amps.append(this_amp)

        # apply the amplitude, if it's meaningful
        if np.abs(this_amp) > 0.001:
            apply_or_not = input('Apply this correction? (Y/n) >>> ')
            if apply_or_not == 'Y':
                tuned = set_dm_shape(n, this_amp, refsrf=refsrf)
            else:
                print('Not applying based on user input.')
                tuned = DM.getSurf()
        else:
            print('Correction too small, not applying this mode')
            tuned = DM.getSurf()

        # save intermediate map
        np.save(os.path.join(date_dir, 'noll' + str(n) + 'map_' + mode + '_sf' + str(fiber) + '_' + date + '.npy'), tuned)
        time.sleep(0.5)
        
        # turn on tracking before pd reads
        track.start_tracking()
        time.sleep(2)

        # Just use ini_bkgd
        # this_bkgd = _take_pd_backgrounds(nreads=50)
        this_pdread = tools.read_pd(nreads).mean() * 1e4 - ini_bkgd
        print('Current PD read: ' + str(this_pdread))

        all_pdread.append(this_pdread)

    final_file = 'NCPA_map_' + mode + '_sf' + str(fiber) + '_' + date + '.npy'
    # np.save(os.path.join(date_dir, final_file), tuned)
    np.save(os.path.join(DM.flatdir, final_file), tuned)

    # Plot the fits, and pd read over time
    for n in modes_to_scan:
        _ = plot_zernike_scan_peaks(n, show=False, date_dir=date_dir)

    startTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
    tstamp = startTime.replace('-','')
    tstamp = tstamp.replace(':','')
    tstamp = tstamp.replace(' ','')
    tstamp = tstamp[2:-2]

    final_pdread = tools.read_pd(nreads).mean() * 1e4
    all_pdread.append(final_pdread)
    print('Final PD read: ' + str(final_pdread))

    fig, ax = plt.subplots()
    str_modes = [str(j) for j in modes_to_scan]
    str_modes = np.insert(str_modes, 0, 'Ini')
    str_modes = np.insert(str_modes, -1, 'Final')

    ax.scatter(str_modes, all_pdread)
    ax.set_xlabel('After noll index')
    ax.set_ylabel('PD power * 1e4 (V)')
    plt.savefig(os.path.join(date_dir, tstamp+'_pdreads.png'), bbox_inches='tight')
    # time.sleep(3)
    # plt.close()
    print('The final DM map is saved as ' + final_file)

    return tuned

def ncpa_scan_tests(start_noll=4, final_noll=12, start=-0.3, stop=0.3, step=0.05, nreads=10, apply_all=True, apply_end=False, date = '221006', refsrf='current', fiber=2,
                grid_size=1, stepsize=2, verbose=False, mode='ds'):

    # Start with previous flat map - try to find automatically. Assume they are saved in ~/dev/
   # all_flats = glob.glob('NCPA_map_' + mode + '*')
    all_flats = glob.glob(os.path.join(DM.flatdir, 'NCPA_map_' + mode + '*'))  
    guess_latest = np.sort(all_flats)[-1]

    print('Suggesting the initial flat map: ' + guess_latest)
    ini_flat_good = input('Is the initial flat map good? (Y/n) >>> ')
    if ini_flat_good == 'Y':
        ini_flat = guess_latest
    else:
        user_ini_flat = input('Please input name of initial flat map >>> ')
        ini_flat = user_ini_flat.strip()
    print('Applying DM Map ' + ini_flat)

    old_flat = np.load(os.path.join(DM.flatdir, ini_flat))
    DM.setSurf(old_flat)

    # initial pd read of background
    ini_bkgd = _take_pd_backgrounds(nreads=50)
    time.sleep(1)

    Path = '/nfiudata/FAM_Scans/'

    #Path = '/home/nfiudev/dev/FAM_Scans/'
    date_dir = os.path.join(Path, date)
    if not os.path.exists(date_dir):
        os.makedirs(date_dir)

    modes_to_scan = np.linspace(start_noll, final_noll, final_noll-start_noll+1, dtype=int)
    all_amps = []
    all_pdread = []
    amps_applied = []

    # start tracking if not already
    if not track.is_tracking():
        track.start_tracking()
        time.sleep(1)
        
    # go to SF2 first
    _acquire_fiber(fiber)
    SAM.set_pos('pd_sf' + str(fiber), block=True)
    time.sleep(2)

    ini_pdread = tools.read_pd(nreads).mean() * 1e4 - ini_bkgd
    all_pdread.append(ini_pdread)
    print('Initial PD read: ' + str(ini_pdread))
    time.sleep(0.5)

    # start scanning modes requested
    for n in modes_to_scan:
        # stop tracking immediately before each scan
        track.stop_tracking()
        this_amp = scan_mode(n, start=start, stop=stop, step=step, grid_size=grid_size, stepsize=stepsize, refsrf=refsrf, fiber=fiber, bkgd=ini_bkgd,
                nreads=nreads, show=True, date_dir=date_dir, verbose=verbose)
        all_amps.append(this_amp)
           
        # Apply every correction w/o prompting, if flag is True
        if apply_all:
            print('Applying correction...')
            tuned = set_dm_shape(n, this_amp, refsrf=refsrf)
        else:
            if apply_end:
            # Wait until the end to apply all correction, if flag is True
                # apply the amplitude, if it's meaningful
                if np.abs(this_amp) > 0.001:
                    apply_or_not = input('Apply this correction? (Y/n) >>> ')
                    if apply_or_not == 'Y':
                        amps_applied.append(this_amp)
                        #tuned = set_dm_shape(n, this_amp, refsrf=refsrf)
                    else:
                        print('Not applying based on user input.')
                        tuned = DM.getSurf()
                else:
                    print('Correction too small, not applying this mode')
                    tuned = DM.getSurf()
            else:
                # apply the amplitude, if it's meaningful
                if np.abs(this_amp) > 0.001:
                    apply_or_not = input('Apply this correction? (Y/n) >>> ')
                    if apply_or_not == 'Y':
                        tuned = set_dm_shape(n, this_amp, refsrf=refsrf)
                    else:
                        print('Not applying based on user input.')
                        tuned = DM.getSurf()
                else:
                    print('Correction too small, not applying this mode')
                    tuned = DM.getSurf()

        if not apply_end:
            # save intermediate map
            np.save(os.path.join(date_dir, 'noll' + str(n) + 'map_' + mode + '_sf' + str(fiber) + '_' + date + '.npy'), tuned)
            time.sleep(0.5)
        
            # turn on tracking before pd reads
            track.start_tracking()
            time.sleep(2)

            # Just use ini_bkgd
            # this_bkgd = _take_pd_backgrounds(nreads=50)
            this_pdread = tools.read_pd(nreads).mean() * 1e4 - ini_bkgd
            print('Current PD read: ' + str(this_pdread))

            all_pdread.append(this_pdread)

    # Apply amps to DM after scan, if flag is True
    if apply_end:
        for n in modes_to_scan:
            for amp in amps_applied:
                tuned = set_dm_shape(n, amp, refsrf=refsrf)

    final_file = 'NCPA_map_' + mode + '_sf' + str(fiber) + '_' + date + '.npy'
    # np.save(os.path.join(date_dir, final_file), tuned)
    np.save(os.path.join(DM.flatdir, final_file), tuned)

    if not apply_end:
        # Plot the fits, and pd read over time
        for n in modes_to_scan:
            _ = plot_zernike_scan_peaks(n, show=False, date_dir=date_dir)

    startTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
    tstamp = startTime.replace('-','')
    tstamp = tstamp.replace(':','')
    tstamp = tstamp.replace(' ','')
    tstamp = tstamp[2:-2]

    if not apply_end:
        final_pdread = tools.read_pd(nreads).mean() * 1e4
        all_pdread.append(final_pdread)
        print('Final PD read: ' + str(final_pdread))

        fig, ax = plt.subplots()
        str_modes = [str(j) for j in modes_to_scan]
        str_modes = np.insert(str_modes, 0, 'Ini')
        str_modes = np.insert(str_modes, -1, 'Final')

        ax.scatter(str_modes, all_pdread)
        ax.set_xlabel('After noll index')
        ax.set_ylabel('PD power * 1e4 (V)')
        plt.savefig(os.path.join(date_dir, tstamp+'_pdreads.png'), bbox_inches='tight')
    # time.sleep(3)
    # plt.close()
    print('The final DM map is saved as ' + final_file)

    return tuned

def sam_scan(sam_y_step = 0.005, scan_range = 0.03, spec_nd = 'nd1', take_spec_exp = False, scan_fib=2):
    # default for spec is nd1
    Light_Src_cmds.set(spec_nd)
    print('Inserted ' + spec_nd)
    
    # load presets 
    SAM.load_presets()
    time.sleep(0.5)
    
    # mmove to on slit
    SAM.set_pos("on_slit", block=True)
    time.sleep(2)

    set_nspec_reads(itime)
    
    # Take spec exposures with fibers 1, 2, 3, 4
    # see if orders are all in the detector. If some are cutoff, move SAM in x position
    if take_spec_exp:
        fibers=[1,2,3,4]
        print('Taking spec exposures to check if all spectral orders are in NIRSPEC.')
        track.start_tracking()
        sam_x_finalized = False

        while not sam_x_finalized:
            fig, axes = plt.subplots(1,4, figsize=(16,4))
            for i, fib in enumerate(fibers):
                # acquire and track
                _acquire_fiber(fib)
                track.start_tracking()
                time.sleep(1)

                # take 1 exposure on each SF
                spec.TakeImage(Coadds=1, tint=itime)  # default 1.5 sec integrations
                this_im = spec.im[0].T[:,::-1]
                # plot these spec images
                axes[i].imshow(this_im, vmin=10, vmax=1000)
                #axes[i].set_title('SF'+str(fib))
            
            fig.suptitle('NIRSPEC frames on SF 1, 2, 3, 4 (left to right)')
            # plt.show()
    
            print('Make sure orders are not cut off from spec, especially in fibers 1 and 4.')
            print('If they are, move SAM in x (increase x to move orders up, decrease x to move down.')
            move_sam_x = input('Move SAM x position? (Y/n) >>> ')

            if move_sam_x == 'Y':
                _sam_x_move = input('Move SAM x by how much? (e.g. -0.02) >>> ')
                sam_delta_x = float(_sam_x_move.strip())

                cur_sam_pos = SAM.get_pos()
                cur_y = cur_sam_pos[1]
                goal_x = cur_sam_pos[0] + sam_delta_x
                SAM.set_pos([goal_x, cur_y], block=True, move=[1,0])
            
            # if we're all good, quit the while loop
            else:
                sam_x_finalized = True
    
    plt.close()
    # start the sam scan in y direction
    ini_y = np.round(SAM.get_pos()[1], 2)
    ini_x = SAM.get_pos()[0]
    print('Initial x y positions:')
    print(ini_x, ini_y)

    start_y = ini_y - scan_range
    stop_y = ini_y + scan_range

    # default is to scan for SF2 only
    print('Scanning SAM in y direction for SF2...')
    best_y = slit_scan(start=start_y, stop=stop_y, step=sam_y_step, fibers=[scan_fib], axis='y')

    # go to optimal position
    time.sleep(2)
    print('Moving to optimal y position.')
    SAM.set_pos([ini_x, best_y], block=True, move = [0,1])

    print('Remember to close the figure after viewing it.')

# evaluates calibrations with NIRSPEC exposures
def evaluate_calib(navg = 5, plot=True, spec_nd='nd1'):
    fibers=[1,2,3,4]
    print('Taking spec exposures on fibers 1, 2, 3, 4.')
    print('Using ND1 and tint 1.5 sec.')

    # default for spec is nd1
    Light_Src_cmds.set(spec_nd)
    print('Inserted ' + spec_nd)

    track.start_tracking()

    all_data = {'1':[], '2':[], '3':[], '4':[]}

    for fib in fibers:
        # acquire and track
        _acquire_fiber(fib)
        track.start_tracking()
        time.sleep(1)

        for i in range(navg):
            # take 1 exposure on each SF
            spec.TakeImage(Coadds=1, tint=itime)  # default 1.5 sec integrations
            this_im = spec.im[0].T[:,::-1]
            all_data[str(fib)].append(this_im)

    bkgd, _, badpixmap = _take_backgrounds(itime, 1, nframes=5)
    which_Order = 2

    mean_flux = []
    for fib in fibers:
        print('Fiber ' + str(fib))
        all_frames = all_data[str(fib)]
        datacube = np.asarray(all_frames)
        # print(datacube.shape)
        trace_locs = _calc_trace(datacube, badpixmap, N_order, N_fiber, plot=plot, bkgd=bkgd)
        # print('trace locs')
        # print(trace_locs.shape, trace_locs)

        fluxes = []
        res_list = parallel_analysis(len(all_frames), all_frames, which_Order, trace_locs, N_order, N_fiber)
        for j in range(len(all_frames)):
            # get order_flux, which is 2nd output (index=1)
            order_flux = res_list[j].get()[1]
            int_flux = res_list[j].get()[0][0]
            # print(res_list[j].get())
            # print('int flux')
            # print(int_flux)
            fluxes.append(int_flux)

        mean_flux.append(np.nanmean(fluxes))
        print(fluxes)

    Path = Organization.get_path('Spec')

    fig, ax = plt.subplots(figsize=(8,4))
    fs = 15
    ax.scatter(fibers, mean_flux, color='black')
    ax.set_xlabel('Science Fiber #',fontsize=fs)
    ax.tick_params(labelsize=fs)
    ax.set_yscale('linear')
    ax.set_ylabel('Integrated Flux on NIRSPEC',fontsize=fs)
    plt.savefig(Path[0] + 'evaluate_calib.png',dpi=200)
    plt.show()

def parallel_analysis(nproc, frames, which_Order, trace_locs, N_order, N_fiber, bkgd=np.zeros((2048,2048)), bpmap=np.ones((2048,2048))):
    # TODO::: cleanup this open pool
    pool = mp.Pool(nproc)
    res_list = []
#    with mp.Pool(nproc) as pool:   
    for fr in frames:
        res_list.append(pool.apply_async(_total_flux_frame, args=(fr, bpmap, N_order, N_fiber,
                    bkgd, trace_locs, which_Order)))

    return res_list
    
def vfn_take_bright_trace(fibers, N_images=3, refsrf='current', date=None, verbose=False):
    '''Do trace finding (assumes a bright source is already set up)
    
    Inputs:
        __(Required)__
        fibers =    list of fibers to use
        __(Optional)__
        N_images =  Number of images to use as samples on each fiber (default 3)
        refsrf =    str: DM surface to use; default 'current'=current DM shape
        date =      str: name of folder for data saving; usually yymmdd (default today's date)
        verbos =    flag whether to print status updates (default False)
    
    Example:
        spec_ims, fluxes_all, trace_locs_all = vfn_take_bright_trace([3,4], verbose=True)
        -- can unpack traces as trace_sf3, trace_sf4 = trace_locs_all
    ''' 
    # Set/Create path for saving results
    Path = '/home/nfiudev/dev/VFN_Cals/'
    if date is None:
        date = datetime.date.today().strftime("%y%m%d")
    date_dir = os.path.join(Path, date)
    if not os.path.exists(date_dir):
        os.makedirs(date_dir)
    
    # Get the right flatmap based on input
    if refsrf != 'current':
        # non-'current' string provide so get and set the desired flat
        flat = DM.setFlatSurf(refsrf)
        if verbose:
            print("'%s' DM map loaded"%refsrf)
    else:
        if verbose:
            print("current DM map used")

    if verbose:
        print('taking backgrounds')
   # Get/Take backgrounds on nirspec
    bkgd, _, bpm = _take_backgrounds(itime, coadds, nframes=5)
    # NOTE: _take_backgrounds turns on and leaves tracking on when it is done...
    if verbose:
        print('backgrounds complete, now taking images for trace finding')
        
   # start tracking if not already (in case they edit _take_backgrounds to not auto-set tracking on)
    if not track.is_tracking():
        track.start_tracking()      
        
    # Ensure input fibers is an iterable (in case a single fiber was provided)
    try: len(fibers)
    except TypeError: fibers = [fibers]
    
    # pre-allocate output variables
    spec_ims = np.empty((len(fibers), N_images, 2048, 2048))    # array of spec images (assuming 2048x2048 format)
    trace_locs_all = []
    fluxes_all = []
    # Iterate over fibers taking images on each
    for fibind, fiber in enumerate(fibers):
        # go to given fiber
        _acquire_fiber(fiber)
        time.sleep(2)
        if verbose:
            print('algined to fiber %d, taking %d nirspec images on this fiber'%(fiber, N_images))
        # Take the desired number of images
        for imind in range(N_images):
            spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='focus_scan_step'+str(i))
            spec_ims[fibind,imind] = spec.im[0].T[:,::-1] - bkgd    
    
        if verbose:
            print('imaging on this fiber complete. Starting trace finding')
        # Use flux finder to simultaneously find trace and verify good flux in frame
            # (remember that not providing loaded_Trace will force it to compute a trace and return it as well)
        if fibind > 0:
            plt.figure()    # create new figure for plotting when more than one fiber is used
        fluxes, trace_locs  = compute_spec_flux(spec_ims[fibind], which_Order=6, bpmap=bpm, bkgd=bkgd, plotTrace=True)
        trace_locs_all.append(trace_locs)
        fluxes_all.append(fluxes)
        
        print('fluxes for fiber {}: {}'.format(fiber, fluxes))

    ## Do we need to save scans array somehow or is it fine to make it an output
    return spec_ims, fluxes_all, trace_locs_all
    

def vfn_ncpa_scan(noll, fiber, loaded_Trace, start=-0.3, stop=0.3, step=0.01, which_Order=6, refsrf='current', N_images=3, date=None, verbose=False):
    '''Inputs: 
        ___(mandatory)___
        noll =          noll index to scan
        fiber =         fiber to scan on
        loaded_Trace =  trace to use (MUST BE PROVIDED since VFN scans will naturally be low-flux so trace extraction will not be reliable)
    
        ___(optional)___
        start, stop, step = elements to make amplitude vector for given noll (in BMC units)
        which_Order =   spectral order for which to compute flux (index-0)
        refsrf =        str: stringname of flat map to use (can be 'current' to use current DM surface)
        N_images =      number of samples to take and median at each amplitude (default 3)
        date =          str: name of output output save directory within ~/dev/VFN_Cals/, generally a date in yymmdd format
        verbose =       flag to denote if code should print status updates
        
    Sample use:
        zpts, spec_ims, fluxes, fit_coeffs = vfn_ncpa_scan(4, 3, trace_sf3, start=-0.2, stop=0.2, step=0.05, which_Order=6, refsrf='current', date=None, verbose=True)
    '''
    # Set/Create path for saving results
    Path = '/home/nfiudev/dev/VFN_Cals/'
    if date is None:
        date = datetime.date.today().strftime("%y%m%d")
    date_dir = os.path.join(Path, date)
    if not os.path.exists(date_dir):
        os.makedirs(date_dir)
    
    # Get the right flatmap based on input
    if refsrf != 'current':
        # non-'current' string provide so get and set the desired flat
        flat = DM.setFlatSurf(refsrf)
        if verbose:
            print("'%s' DM map loaded"%refsrf)
    else:
        # 'current' provided so register the current DM surface as the flat
        flat = DM.getSurf()
        if verbose:
            print("current DM map used")

    try:
        if verbose:
            print('taking backgrounds')
       # Get/Take backgrounds on nirspec
        bkgd, _, bpm = _take_backgrounds(itime, coadds, nframes=5)
        # NOTE: _take_backgrounds turns on and leaves tracking on when it is done...
        if verbose:
            print('backgrounds complete, now aligning to fiber')
            
       # start tracking if not already (in case they edit _take_backgrounds to not auto-set tracking on)
        if not track.is_tracking():
            track.start_tracking()        
        # go to given fiber
        _acquire_fiber(fiber)
        time.sleep(2)
        # Now that fiber is acquired, stop tracking
        track.stop_tracking()
        if verbose:
            print('algined to fiber, starting zernike scan')

        # Pre-allocate scan points and data arrays
        zpts = np.arange(start, stop+step/2., step)
        spec_ims = np.empty((zpts.size, N_images, 2048, 2048))    # array of spec images (assuming 2048x2048 format)
        #cred2_ims =  np.empty((zpts.size, *tc.Img.get_data(reform=True).shape)) # cred2 ims using current size
        
        # Iterate through amplitudes
        for i in range(zpts.size):
            if verbose:
                print('Amplitude:', np.round(zpts[i],3))
            ### Apply zernike
            dmmap = DM.pokeZernike(zpts[i], noll, bias=flat)
            # Take desired number of samples
            for j in range(N_images):
                ### Read nirspec
                spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='focus_scan_step'+str(i))
                spec_ims[i,j] = spec.im[0].T[:,::-1] - bkgd    
            ### Set DM back where it started to make sure we have something good to track on
            DM.setSurf(flat)
            ### re-acquire fiber to make sure we are well-aligned
            _acquire_fiber(fiber, toggle_track=True)
            time.sleep(0.5)
        
        if verbose:
            print('zernike scan complete. Computing integrated flux in order %d'%which_Order)
        
        # Take median of samples at each amplitude
        spec_ims_med = np.median(spec_ims, axis=1)
            
        # Extract integrated fluxes
            #*** TODO::: modify _total_flux_frame to also return full spectrum (see vfnserver for how to do this)
            # Would need to do this in a non-invasive way. could be a flag that chooses which return version to use with default for normal return
        fluxes = compute_spec_flux(spec_ims_med, which_Order, loaded_Trace, bpm, bkgd)
        
        # Perform second-order polynomial fit (parabolic)
            # TODO::: consider switching to a gaussian which is the more-likely true form
        coeffs = np.polyfit(zpts,fluxes,2)
        
        # Display results and fit
        plt.scatter(zpts,fluxes)
        plt.plot(zpts,[coeffs[0]*z**2+coeffs[1]*z+coeffs[2] for z in zpts])
        plt.xlabel('Amplitudes [BMC units]')
        plt.ylabel('Integrated flux over order %d [cts]'%which_Order)
        plt.title('Fiber %d - Noll %d'%(fiber, noll))
        
        # Find which amplitude minimized the null
        minIndex = fluxes.index(min(fluxes))
        amplitude = zpts[minIndex]
        
        print('==>> Amp. on Noll %d w/ min flux in order %d for fiber %d: %f'%(noll, which_Order, fiber, amplitude))
    except Exception as e:
        print('Error encountered, exiting at this point')
        print('Error text: {}'.format(e))
        return spec_ims, zpts, None, None        

    ## Do we need to save scans array somehow or is it fine to make it an output
    return zpts, spec_ims, fluxes, coeffs

def compute_spec_flux(steps, which_Order=0, loaded_Trace=np.zeros(0), bpmap=np.ones((2048,2048)), bkgd=np.zeros(0), plotTrace=False):
    '''
    Inputs:
    steps: shape (Nframes x 2048 x 2048) of input datacube. 
    which_Order: which order to report integrated flux for
    loaded_Trace: trace to use. If not given, computes trace from input datacube
    bpmap: badpixmap to use, all ones to no badpix correct (default)
    bkgd: background to subtract, all 0 to not subtract (default)
    plotTrace: whether to plot trace. Could be True for diagnostic purposes


    Returns:
    fluxes: list of integrated fluxes in which_Order for each of the input frames
    trace_locs: (if loaded_Trace is not provided) the trace locations
    '''

    fluxes = []
    all_frames = []

    # reorganize the datacube into a list of frames...
    for frame in steps:
        all_frames.append(frame)

    # calc trace if not given one
    if loaded_Trace.shape[0] == 0:
        # use median of the frames to compute trace.
        datacube = np.asarray(all_frames)
        trace_locs = _calc_trace(datacube, bpmap, N_order, N_fiber, plot=plotTrace, bkgd=bkgd)
    else:
        trace_locs = loaded_Trace

    # res_list = []     # I don't think this is necessary, comment out to test, then remove if possible
    if len(all_frames) > 1:
        res_list = parallel_analysis(len(all_frames), all_frames, which_Order, trace_locs, N_order, N_fiber, bkgd, bpmap)
        for j in range(len(all_frames)):
            # get order_flux, which is 2nd output (index=1)
            order_flux = res_list[j].get()[1]
            fluxes.append(order_flux)
    else:
        this_frame = all_frames[0]
        
        # int_flux is the total flux across all the orders
        # order_flux is integrated flux of only 1 order
        int_flux, order_flux = _total_flux_frame(this_frame, bpmap, N_order, N_fiber, bkgd, trace_locs, which_Order)
        fluxes.append(order_flux)

    print(fluxes)
    if loaded_Trace.shape[0] == 0:
        return fluxes, trace_locs
    else:
        return fluxes
