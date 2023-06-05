import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.pylab as pl
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm, Normalize
import matplotlib as mpl
from astropy.io import fits, ascii
from astropy.modeling import models
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
from Nirspec_cmds import Spec_cmds, Scam_cmds, Light_Src_cmds
from Star_Tracker_cmds import Tracking_cmds
from Track_Cam_cmds import TC_cmds
from FAM_cmds import FAM_cmds
from SAM_cmds import SAM_cmds
from PSM_cmds import PSM_cmds
from Telescope_Simulator_cmds import Telescope_Simulator_cmds
import Acquisition
import ktl
import Organization
from PSF_finder import moments, Gaussian2D, PSF_Finder
import sys
import time
import warnings
import datetime
import pdb
from time import sleep

from scipy.signal import medfilt as medfilt

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

print('Current save directory:', Organization.get_path()[0])

cam   = TC_cmds()
track = Tracking_cmds()
TTM   = FAM_cmds()
spec  = Spec_cmds()
scam  = Scam_cmds()
PSM   = PSM_cmds()
SAM   = SAM_cmds()
# connect to flip mirror keyword
flip  = ktl.Service("kpic")["FEUFLPST"]
flip.monitor()
#fib = FEU_Fiber_cmds()
IS = Telescope_Simulator_cmds()
# Light_Src_cmds.set('nd2') ### I think this is what we use but need to check

### These are some parameters that probably don't need to change often
N_order = 7
N_fiber = 1
nskip = 4 ### Only fit every 4 points for speed when using KPIC DRP, still give good enough trace after smoothing
itime = 1.5*1000
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
        time.sleep(0.1)

    return np.mean(A)

def _acquire_fiber(fibr):
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
            sys.stdout.write(u"\u001b[2K")
            sys.stdout.write("\rfiber {} | acquiring".format(fibr)+dots)
            sys.stdout.flush()
            # get location of psf
            params = track.get_psf_parameters()
            # get target location
            goal = track.get_psf_goal()

            # if we're close enough to goal, break
            if params[0] and abs(params[2] - goal[1]) < .5 and abs(params[3] - goal[2]) < .5:
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
        bkg_names.append(Path[0]+'bkgd_'+str(i)+'.fits')
        if not os.path.isfile(bkg_names[-1]):
            retake = True

    if retake:
        ### Stop tracking and move off-slit for backgrounds
        if track.is_tracking():
                track.stop_tracking()
                # Wait for the Tracking loop to complete its last iteration
                time.sleep(2)
        ### Get initial TTM position
        TTM_CP = TTM.get_pos()
        ### Move TTM and wait for new position
        TTM.set_pos([9000, 9000], block=True)

        ### Take backgrounds
        for i in range(nframes):
            ### Take backgrounds, save for DRP 
            spec.TakeImage(Coadds=coadds, tint=itime, Save=True, filename='bkgd_'+str(i))
        TTM.set_pos([TTM_CP[0], TTM_CP[1]], block=True)

    ### Return to original TTM position and restart tracking
    bkgd, smoothed_thermal_noise, badpixmap = background.make_badpixmap(bkg_names)

    track.start_tracking()

    return bkgd, smoothed_thermal_noise, badpixmap

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

def FEU_Scan_XY(xstart=None, xstop=None, xstep=0.4, ystart=None, ystop=None, ystep=0.4, Scam_tint=3000, usavename='', plot=True):
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
    Track_Flag = track.is_tracking()
    goal0 = track.get_goal()[0]
    try:
        SAM.set_pos("off_slit")

        # Prepare scanned position
        List_pos_x = np.arange(xstart, xstop+xstep/2., xstep)
        List_pos_y = np.arange(ystart, ystop+ystep/2., ystep)

        # print(List_pos_x, List_pos_y)

        # Compute the total number of steps
        nbsteps = np.size(List_pos_x)*np.size(List_pos_y)

        if Track_Flag:
            track.stop_tracking()
            # Wait for the Tracking loop to complet its last iteration
            time.sleep(20 * cam.get_tint())
            print('Tracking loop gain set to 0 for the scan.')
        
        # Get the current position of the TTM
        TTM_CP = TTM.get_pos()

        # Move the TTM 
        TTM.set_pos([9000,9000], block = True)
        # Take Scam calibration images
        scam.TakeImage(tint=Scam_tint)
        scam_bkgd = scam.im[0]
                
        # Return to initial TTM position
        TTM.set_pos([TTM_CP[0],TTM_CP[1]], block=True)
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
                tmp_im = medfilt(scam.im[0]-scam_bkgd,3)
                x,y = np.where(tmp_im == np.max(tmp_im))
                x = int(x[0])
                y = int(y[0])
                ### Get counts in area around max flux on SCAM
                flux = np.sum(tmp_im[x-3:x+4,y-3:y+4])
                Results[0,tmp_x,tmp_y] = flux #/cred2flux
                # Position of the Zaber X
                Results[1,tmp_x,tmp_y] = pos_x
                # Position of the Zaber Y
                Results[2,tmp_x,tmp_y] = pos_y  

                # Compute average time per iteration
                avr_time = (time.time()-Time_ini)/(ite_nb+1.)
                # Compute time before end of scan
                timeleft = round((nbsteps-ite_nb)*avr_time,0)
                # Increment ite number
                ite_nb += 1
                # Prepare information to print
                text  = 'Remaining Time %05.d sec' %(timeleft)
                text += ' -- %03d/%03d ite left' %(nbsteps-ite_nb,nbsteps)
                text += ' -- Flux = %06.3f' %(Results[0,tmp_x,tmp_y])
                # Print the information in the terminal
                sys.stdout.write('\r FEU_Scan_XY: ' + text)
                sys.stdout.flush()

        # Return to initial TTM  and Zaber positions
        TTM.set_pos([TTM_CP[0],TTM_CP[1]], block=True)
        SAM.set_pos([SAMX_ini, SAMY_ini], block=True)
        PSM.set_pos([PSMX_ini, PSMY_ini], block=True)
        if Track_Flag: track.start_tracking()
        
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
            fig = _plot_fiu_ttm_scan(Rx, Ry, Flux.T, fitmod.T, residual.T, offset_x, offset_y, (Path,filename), usavename)
            fig.show()

        return offset_x, offset_y

    except (Exception, KeyboardInterrupt) as e:
        print('Exception encountered, resetting initial state')
        TTM.set_pos([TTM_CP[0],TTM_CP[1]], block=True)
        SAM.set_pos([SAMX_ini, SAMY_ini], block=True)
        PSM.set_pos([PSMX_ini, PSMY_ini], block=True)
        if Track_Flag: track.start_tracking()
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
    F_ini = PSM.get_foc_pos()[0]
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

        datacubes = []
        for fibr in fibers:
            print('Fiber', fibr)
            # add extra line for \033[F\033[K to delete
            print()
            _ = _acquire_fiber(fibr)
            foci = np.arange(start, stop+step/2., step)
            datacube = np.empty((foci.size, 2048, 2048))
            ipos = PSM.get_foc_pos()[0]
            for i in range(foci.size):
                try:
                    ### Move FEU
                    PSM.set_foc_pos(foci[i], block=True)

                    ### Get and write current position to terminal
                    tmp = PSM.get_foc_pos(update = False)[0]
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
                PSM.set_foc_pos(F_ini, block=True)
            except:
                print('PSM focus failure, position: ', ipos)
            datacubes.append(datacube)

        if len(datacubes)==1:
            return foci, datacubes[0]
        else:
            return foci, datacubes

    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        PSM.set_foc_pos(F_ini, block=True)
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
            plt.ylabel('Flux [DN]')
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

def _take_slit_scan(start=3.275, stop=3.325, step=0.005, bkgd=None, Path=None, savename='', fibers=None, save=True):
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
                SAM.set_pos([xpts[i], ipos[1]], block=True, move=[1,0])

                tmp = SAM.get_pos(update = False)
                txt = 'X = %5.2f - Y = %5.2f' %(tmp[0], tmp[1])
                sys.stdout.write('\rFEU TTM Position: ' + txt)
                sys.stdout.flush()

                ### Take an image, get data-background
                spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='slit_scan_step'+str(i))
                datacube[i] = spec.im[0].T[:,::-1] - bkgd
            
            if save:
                usavename = '_fiber'+str(fibr)+savename
                hdu = fits.PrimaryHDU(datacube)
                xstr = '%5.2f | '*xpts.size
                hdu.header['xpts'] = xstr%tuple(xpts)
                hdu.header['fiber'] = fibr
                hdu.writeto(Path[0] + Path[1] + '_slit_scan_cube'+usavename+'.fits', overwrite=True)
            print('')
            ## Set FEU back to initial position
            SAM.set_pos(ipos, block=True)
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

def _analyze_slit_scan(datacube=None, xpts=None, badpixmap=np.ones((2048,2048)), plot=True, Path=None, fibers=None, savename='', plotTrace=False):
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
        int_flux = _total_flux_fast(datacube, badpixmap, N_order, N_fiber, plotTrace=plotTrace)
        transmissions*=int_flux/np.nanmax(int_flux)

        ### Output and plot
        best_pos = np.round(xpts[np.argmax(int_flux)],3)
        print('Best FEU TTM x position:', best_pos, ustr)
        np.savez(Path[0]+Path[1]+'_slit_scan_data'+usavename+'.npz', int_flux=int_flux)
        if plot:
                plt.plot(xpts, int_flux, label='Fiber '+str(fibr))  
    if plot:
        plt.axvline(x=xpts[np.argmax(transmissions)], color='r', linestyle='--',label='Best FEU TTM Position')
        plt.text(xpts[np.argmax(transmissions)]+0.005, 0.75*np.max(int_flux), str(np.round(xpts[np.argmax(transmissions)],3)))
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
    if len(fibers) > 1:
        print('FEU Position for maximum transmission in all fibers:', np.round(xpts[np.argmax(transmissions)],3))
        return np.round(xpts[np.argmax(transmissions)],3), np.nanmax(transmissions)
    else:
        return np.round(xpts[np.argmax(int_flux)],3)

def slit_scan(start=3.275, stop=3.325, step=0.005, plot=True, fibers=None, plotTrace=False):
    '''
    Scan over given FEU TTM region to maximize flux
    '''
    Path = Organization.get_path('Spec')
    Path = list(Path)
    Path[1] =Path[1].split('.')[0]
    ### Take backgrounds
    bkgd, smoothed_thermal_noise, badpixmap = _take_backgrounds(itime, coadds, nframes=3)
    badpixmap = np.ones((2048,2048))
    # badpixcube = np.expand_dims(badpixmap, axis=0)
    # bkgd = np.zeros((2048, 2048))
    # smoothed_thermal_noise = np.zeros((2048, 2048))

    ### Set up scanning points
    xpts, datacube = _take_slit_scan(start=start, stop=stop, step=step, bkgd=bkgd, Path=Path, fibers=fibers)

    best_pos = _analyze_slit_scan(datacube=datacube, xpts=xpts, badpixmap=badpixmap, plot=plot, Path=Path, fibers=fibers, plotTrace=plotTrace)

    return best_pos

def _take_fiu_ttm_scan(start=-4, stop=4, step=1, fibers=None, bkgd=None, Path=None, savename=''):
    '''
    Takes 2D FIU TTM scan
    '''
    Track_Flag = track.is_tracking()
    if Track_Flag:
        goal0 = track.get_goal()[0]
    try:
        if bkgd is None:
            bkgd, _, _, = _take_backgrounds(itime, coadds, nframes=5)
        if fibers is None:
            fibers = [goal0]
        if Path is None:
            Path = Organization.get_path()
        xpts = np.arange(start, stop+step/2., step)
        ypts = np.arange(start, stop+step/2., step)
        datacubes = []
        cred2fluxes = []
        for fibr in fibers:
            cred2fluxes = []
            print('Fiber', fibr)
            usavename = '_fiber'+str(fibr)+savename
            _ = _acquire_fiber(fibr)
            datacube = np.empty((xpts.size*ypts.size, 2048, 2048))
            ### Scan over 2D array
            cnt = 0
            nframes = str(int(xpts.size*ypts.size))
            t0 = datetime.datetime.now()
            for i in range(xpts.size):
                for j in range(ypts.size):
                    ### Move FEU TTM and wait to update
                    track.set_user_offset(offset_x=xpts[i], offset_y=ypts[j])
                    track.wait_on_target()

                    t1 = datetime.datetime.now()
                    # perit = float(t1-t0)/float(cnt+1)
                    # tleft = np.round((xpts.size*ypts.size-cnt+1)*perit,3)
                    tmp = track.get_psf_goal()
                    # txt = 'X = %5.2f - Y = %5.2f, Time remaining: %5.2f s, Frame %3.0f / '%(float(tmp[1]), float(tmp[2]), tleft, float(cnt+1)) + nframes
                    txt = 'X = %5.2f - Y = %5.2f, Frame %3.0f / '%(float(tmp[1]), float(tmp[2]), float(cnt+1)) + nframes
                    sys.stdout.write('\rFIU Pixel Targets: ' + txt)
                    sys.stdout.flush()

                    ### Take an image, get data-background
                    spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='slit_scan_step'+str(i))
                    cred2fluxes.append(_get_cred2_flux())
                    datacube[cnt] = spec.im[0].T[:,::-1] - bkgd
                    cnt +=1

            cred2fluxes = np.asarray(cred2fluxes)
            hdu = fits.PrimaryHDU(datacube)
            xstr = '%5.2f | '*xpts.size
            ystr = '%5.2f | '*ypts.size
            cstr = '%5.2f | '*cred2fluxes.size
            hdu.header['fiber'] = fibr
            hdu.header['xpts'] = xstr%tuple(xpts)
            hdu.header['ypts'] = ystr%tuple(ypts)
            hdu.header['cred2flx'] = cstr%tuple(cred2fluxes)
            hdu.writeto(Path[0] + Path[1] + '_fiu_scan_cube'+usavename+'.fits', overwrite=True)
            print('')
            ## Set FEU back to initial position
            track.set_user_offset(offset_x=0, offset_y=0)
            datacubes.append(datacube)

        if len(datacubes) == 1:
            return xpts, ypts, cred2fluxes, datacubes[0]
        else:
            return xpts, ypts, cred2fluxes, datacubes

    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        track.set_user_offset(offset_x=0, offset_y=0)
        if Track_Flag: 
            track.set_goal(goal0)
            track.start_tracking()
        print('Reset complete')
        # print(e)
        raise e
        assert 1==0
        return np.asarray([np.nan]), np.asarray([np.nan]), np.asarray([np.nan]), np.nan*np.zeros((2048,2048))

def _analyze_fiu_ttm_scan(datacube=None, xpts=None, ypts=None, cred2fluxes=None, badpixmap=np.ones((2048,2048)), plot=True, Path=None, fibers=None, savename='', fit=True, overwrite=False, plotTrace=False):
    if Path is None:
        Path = ('', '')
    if fibers is None:
        fibers = [track.get_goal()[0]]
    ### Use a list of cubes for multiple fibers
    crdflxs = cred2fluxes
    if datacube is not None:
        datacube = np.asarray(datacube)
    datacubes = []
    if datacube is None:
        for i in range(len(fibers)):
            cname = input('Enter datacube name, fiber '+str(int(fibers[i]))+' >> ').replace('\'', '').replace(' ', '')
            fname = Path[0]+Path[1]+cname
            if not os.path.exists(fname):
                print('Invalid datacube filename!')
                print(fname)
                return -1
            else:
                datacubehdu = fits.open(fname)
                xpts = np.asarray(datacubehdu[0].header['xpts'].split('|')[:-1],dtype=float)
                ypts = np.asarray(datacubehdu[0].header['ypts'].split('|')[:-1],dtype=float)
                # crdflxs = np.asarray(datacubehdu[0].header['cred2flx'].split('|')[:-1],dtype=float)
                # print(xpts)
                # xpts = np.arange(-4,4.5,1)
                # ypts = np.arange(-4,4.5,1)
                crdflxs = np.asarray(datacubehdu[0].header['cred2flx'].split('|')[:-1],dtype=float)
                datacubes.append(datacubehdu[0].data)
                Path = ('', fname.split('_fiu')[0])
                # print(Path)
    elif datacube.ndim==3:
        datacubes.append(datacube)
    else:
        datacubes = datacube
    datacubes = np.asarray(datacubes)
    ### pts are steps if not given
    if xpts is None or ypts is None:
        xpts = np.arange(int(np.sqrt(len(datacube))))
        ypts = np.arange(int(np.sqrt(len(datacube))))
        ustr = 'step'
    else: ustr = 'pos'
    if crdflxs is None:
        crdflxs = np.ones(len(datacube))
        # crdflxs = np.asarray(datacubehdu[0].header['cred2flx'].split('|')[:-1],dtype=float)
        # print(crdflxs)
    # print(crdflxs)

    bestxs, bestys = [], []
    for i, datacube in enumerate(datacubes):
        usavename = '_fiber'+str(int(fibers[i]))+savename
        print('Fiber', str(int(fibers[i])))
        
        flux_lst = _total_flux_fast(datacube, badpixmap, N_order, N_fiber, plotTrace=plotTrace)
        ### Get total flux in each frame
        int_flux = np.empty((ypts.size, xpts.size))
        crd_flux = np.empty((ypts.size, xpts.size))
        cnt = 0
        for i in range(xpts.size):
            for j in range(ypts.size):
                int_flux[i,j] = flux_lst[cnt]
                crd_flux[i,j] = crdflxs[cnt]
                cnt+=1
        int_flux = int_flux[:,::-1] ### Flip ordering so lower values are at bottom of image
        crd_flux = crd_flux[:,::-1]
        # print(int_flux)
        int_flux/=crd_flux
        # print(int_flux) 
        np.savez(Path[0]+Path[1]+'_fiu_ttm'+usavename+'_int_fluxes.npz', int_flux=int_flux, xpts=xpts, ypts=ypts, crd_flux=crd_flux)

        besty, bestx = np.unravel_index(np.argmax(int_flux), int_flux.shape)
        besty = ypts[besty] ### Need to flip to match with image
        bestx = xpts[bestx]

        flx_ratio = int_flux/itime/crd_flux
        # plt.imshow(flx_ratio)
        # plt.show()

        if not fit:
            print('Best offsets: X = %3.0f Y = %3.0f'%(bestx,besty))
            if plot:
                extent = [xpts[0]-0.5, xpts[-1]+0.5, ypts[0]-0.5, ypts[-1]+0.5]
                plt.imshow(int_flux/np.max(int_flux), extent=extent)
                plt.colorbar()
                plt.title('FIU TTM Offset Scan')
                plt.xlabel('X pixel offset')
                plt.ylabel('Y pixel offset')
                plt.savefig(Path[0]+Path[1]+'_fiu_ttm_'+usavename+'.png', bbox_inches='tight')
                plt.show()
            bestxs.append(bestx)
            bestys.append(besty)

        else:
            Flux  = int_flux
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
            p,u = opt.curve_fit(Gaussian2D, dim, Vec, p0 = SP)

            fitmod = np.reshape(Gaussian2D(dim,p[0],p[1],p[2],p[3],p[4],p[5],p[6]),dim) 
            residual = Flux - fitmod
            #
            offset_x = fx(p[1])
            offset_y = fy(p[2])
            #
            print('Gaussfit:')
            print('Maximum: %8.6f' %(p[0]))
            print('Best position (fixed): X = %06.2f -- Y = %06.2f' %(offset_x,offset_y))
            bestxs.append(offset_x)
            bestys.append(offset_y)

            if plot:
                fig = _plot_fiu_ttm_scan(Ry, Rx, Flux.T, fitmod.T, residual.T, offset_x, offset_y, Path, usavename)
                fig.show()
        ### Update fibers if given a list of fibers
        # if overwrite and fibers[i] !=-1:
        #   if np.abs(bestxs[i]) > 0.3 or np.abs(bestys[i]) > 0.3:
        #       _,_, _, icords = track.get_goal()
        #       xi, yi = icords 
        #       track.set_fiber_loc((xi+bestxs[i], yi+bestys[i]),fibers[i])

    if len(bestxs) == 1:
        return bestxs, bestys
    else:
        return bestxs, bestys

def _plot_fiu_ttm_scan(Rx, Ry, Flux, model, residual, offset_x, offset_y, Path, savename):
    fig = plt.figure(num = 1, figsize = (15,4))
    
    gs1 = GridSpec(1, 6)
    gs1.update(left=0.05, right=0.95, top = 0.84, bottom = 0.05, hspace=0.4, wspace=0.4)

    # Create a subplot for the image
    plt.subplot(gs1[:,:2])
    plt.set_cmap('RdBu')
    # Add the title to this sub image
    plt.title('Injection Map')
    # Plot limits
    limits = [np.min(Rx)-0.5,np.max(Rx)+0.5,np.min(Ry)-0.5,np.max(Ry)+0.5]#FEU_Sc
    # limits = [np.min(Rx),np.max(Rx),np.min(Ry),np.max(Ry)]
    # The injection map is display in linear scale
    plt.imshow(np.abs(Flux),extent=limits, vmin=-1, vmax=1)
    # Modify the axis: one ticks every 3 pixel
    # X_ticks = np.arange(limits[0]+0.5, limits[1], 1)
    X_ticks = np.arange(limits[0], limits[1]+1, 1)
    plt.gca().set_xticks(X_ticks)
    # Y_ticks = np.arange(limits[2]+0.5, limits[3], 1)
    Y_ticks = np.arange(limits[2], limits[3]+1, 1)
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
    X_ticks = np.arange(limits[0], limits[1]+1, 1)
    plt.gca().set_xticks(X_ticks)
    # Y_ticks = np.arange(limits[2]+0.5, limits[3], 1)
    Y_ticks = np.arange(limits[2], limits[3]+1, 1)
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
    X_ticks = np.arange(limits[0], limits[1]+1, 1)
    plt.gca().set_xticks(X_ticks)
    # Y_ticks = np.arange(limits[2]+0.5, limits[3], 1)
    Y_ticks = np.arange(limits[2], limits[3]+1, 1)
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
        
    plt.savefig(fullname, bbox_inches='tight', pad_inches=0.25, dpi=600)
    return fig
    # plt.show()

def fiu_ttm_scan(start=-4, stop=4, step=1, fibers=None, savename='', plot=True, fit=True, overwrite=False, plotTrace=False):
    '''
    Do FIU_TTM scan, start and stop tell how many pixels to offset and step sets bounds
    Scan is 2D, plots a map of integrated flux vs offset
    '''
    Path = Organization.get_path('Spec')
    Path = list(Path)
    Path[1] = Path[1].split('.')[0]
    ### Take backgrounds
    bkgd, smoothed_thermal_noise, badpixmap = _take_backgrounds(itime, coadds, nframes=3)
    badpixmap = np.ones((2048,2048))
    # badpixcube = np.expand_dims(badpixmap, axis=0)
    # bkgd = np.zeros((2048, 2048))
    # smoothed_thermal_noise = np.zeros((2048, 2048))

    xpts, ypts, cred2fluxes, datacube = _take_fiu_ttm_scan(start=start, stop=stop, step=step, fibers=fibers, bkgd=bkgd, Path=Path, savename=savename)

    xoffset, yoffset = _analyze_fiu_ttm_scan(datacube=datacube, xpts=xpts, ypts=ypts, cred2fluxes=cred2fluxes, badpixmap=badpixmap, plot=plot, fit=fit, Path=Path, plotTrace=plotTrace, savename=savename, fibers=fibers, overwrite=overwrite)

    return xoffset, yoffset

def _take_ktl_scan(device='nsmotor', param='rotatorval', stop=-2.0, start=0., step=0.25, bkgd=None, Path=None, savename=''):
    '''
    Only does the currently tracked fiber, if we wind up needing it more I'll update the code
    '''
    ipos = ktl.read(device, param)
    Track_Flag = track.is_tracking()
    if Track_Flag:
        goal0 = track.get_goal()[0]
    try:
        if Path is None:
            Path = Organization.get_path('Spec')
            Path = list(Path)
            Path[1] =Path[1].split('.')[0]
        if bkgd is None:
            bkgd, _, _ = _take_backgrounds(itime, coadds, nframes=3)

        xpts = np.arange(start, stop+step/2, step)
        datacube = np.empty((rotpts.size, 2048, 2048))
        ipos = ktl.read(device, param)
        for i in range(len(xpts)):
            ### Move the rotator
            ktl.write(device, param, xpts[i])
            time.sleep(2) ### Not sure if blocking

            tmp = ktl.read(device, param)
            txt = '%5.2f'%tmp
            sys.stdout.write('\r'+param+': ' + txt)
            sys.stdout.flush()

            ### Take an image, get data-background
            spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='slit_scan_step'+str(i))
            datacube[i] = spec.im[0].T[:,::-1] - bkgd
        ### Reset rotator and save cube
        ktl.write(device, param, ipos)
        hdu = fits.PrimaryHDU(datacube)
        hdu.writeto(Path[0] + Path[1] + '_'+param+'_scan_cube'+savename+'.fits', overwrite=True)
        print('')

        return xpts, datacube

    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        ktl.write(device, param, ipos)
        if Track_Flag: 
            track.set_goal(goal0)
            track.start_tracking()
        print('Reset complete')
        # print(e)
        raise e
        return np.asarray([np.nan]), np.nan*np.zeros((2048,2048))

def _analyze_ktl_scan(datacube=None, xpts=None, badpixmap=np.ones((2048,2048)), plot=True, Path=None, plotTrace=False):
    ### Check inputs
    if Path is None:
        Path = ('', '')
    if datacube is None:
        cname = input('Enter datacube name >> ').replace('\'', '').replace(' ', '')
        fname = Path[0]+Path[1]+cname
        if not os.path.exists(fname):
            print('Invalid datacube filename!')
            print(fname)
            return -1
        else:
            datacube = fits.open(fname)[0].data
            Path = ('', fname.split('_'+param)[0])
    if xpts is None:
        xpts = np.arange(len(datacube))

    ### Get total flux
    int_flux = _total_flux_fast(datacube, badpixmap, N_order, N_fiber, plotTrace=plotTrace)

    print('Best '+param+' value:', np.round(rotpts[np.argmax(int_flux)],3))
    np.savez(Path[0]+Path[1]+'_'+param+'_scan_data.npz', int_flux=int_flux)
    if plot:
        plt.plot(xpts, int_flux)
        plt.ylabel('Total Flux on Spec')
        plt.xlabel(param)
        plt.savefig(Path[0]+Path[1]+'_'+param+'_scan.png', bbox_inches='tight')
        plt.show()

    return np.round(rotpts[np.argmax(int_flux)],3)

def ktl_scan(device='nsmotor', param='rotatorval', start=-2.0, stop=0., step=0.25, plot=True, plotTrace=False):
    '''
    Scan through whatever you want, as long as it's in ktl. Plots integrated flux vs parameter value
    '''
    Path = Organization.get_path('Spec')
    Path = list(Path)
    Path[1] =Path[1].split('.')[0]

    bkgd, smoothed_thermal_noise, badpixmap = _take_backgrounds(itime, coadds, nframes=3)
    badpixmap = np.ones((2048,2048))

    xpts, datacube = _take_ktl_scan(device=device, param=param, start=start, stop=stop, step=step, bkgd=bkgd, Path=Path)

    xbest = _analyze_ktl_scan(datacube=datacube, xpts=xpts, badpixmap=badpixmap, plot=plot, Path=path, plotTrace=plotTrace)

    return xbest

def _take_slit_rot_scan(start=-0.35, stop=-0.25, step=0.01, fibers=[1,2,3,4], rstart=None, rstop=None, rstep=None, savename='', Path=None):
    '''
    Rotator commands are probably going to be rewritten when the libraries are updated, but we need this for the 
    next few weeks
    '''

    assert start < stop 
    if rstart is not None:
        assert rstart < rstop
        assert rstep is not None

    irot = ktl.read(device='nsmotor', param='rotatorval')
    ipos = FEU_TTM.get_pos()
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

        ### Make input structures
        if rstart is not None:
            rot_pts = np.arange(rstart, rstop+rstep/2., rstep)
        else:
            rot_pts = np.asarray([0])

        datacubes = np.empty((len(rot_pts)*len(fibers), 2048, 2048))
        rot_init = ktl.read(device='nsmotor', param='rotatorval')
        for j,fibr in enumerate(fibers):
            ### Switch to the fiber
            print('Fiber', fibr)
            usavename = '_fiber'+str(fbr)+savename
            _ = _acquire_fiber(fibr)

            ### Step through rotator points
            for i in range(len(rot_pts)):
                ktl.write('nsmotor', 'rotatorval', rot_pts[i])
                print('Rotator Position:',np.round(rot_pts[i],3))
                ### At each rotator point, take a slit scan
                xpts, slitarr = _take_slit_scan(start=start, stop=stop, step=step, bkgd=bkgd, save=False)
                datacubes[i*len(fibers):i*len(fibers)+j] = slitarr

        hdu = fits.PrimaryHDU(datacube)
        rotstr = '%5.2f | '*rot_pts.size
        xstr = '%5.2f | '*xpts.size
        hdu.header['rotpts'] = rotstr%tuple(rot_pts)
        hdu.header['xpts'] = xstr%tuple(xpts)
        hdu.writeto(Path[0] + Path[1] + '_slit_rot_scan_cube'+savename+'.fits', overwrite=True)
        print('')

        return xpts, rot_pts, datacubes

    except (Exception, KeyboardInterrupt) as e:
        print('\nUnhandled exception encountered, resetting initial state')
        ktl.write(ipos, device='nsmotor', param='rotatorval')
        FEU_TTM.set_pos(ipos, block=True)
        if Track_Flag: 
            track.set_goal(goal0)
            track.start_tracking()
        print('Reset complete')
        # print(e)
        raise e
        return np.asarray([np.nan]), np.asarray([np.nan]), np.nan*np.zeros((2048,2048))

def _analyze_slit_rot_scan(datacubes=None, start=-0.35, stop=-0.25, step=0.005, fibers=[1,2,3,4], rstart=None, rstop=None, rstep=None, plot=True,savename=''):
    '''
    Reduce scan above. Haven't used enough to justify nice header updates
    '''

    ### Check input validity
    assert start < stop 
    if rstart is not None:
        assert rstart < rstop
        assert rstep is not None
    if bkgd is None:
        bkgd, _, _, = _take_backgrounds(itime, coadds, nframes=5)
    if Path is None:
        Path = Organization.get_path('Spec')
        Path = list(Path)
        Path[1] =Path[1].split('.')[0]
    if fibers is None:
        fibers = [-1]

    ### Get the path, drop the fractional seconds
    ### Local override for testing

    ### Make input structures
    if rstart is not None:
        rot_pts = np.arange(rstart, rstop+rstep/2., rstep)
    else:
        rot_pts = np.asarray([0])
    cubes = []
    xpts = np.arange(start, stop+step/2, step)
    best_ttms = np.empty(rot_pts.size)
    max_trans = np.empty(np.size(rot_pts))

    for i in range(len(rot_pts)):
        if datacubes is not None:
            rotcube = datacubes[i*len(fibers):i*len(fibers)+len(fibers)]
        else:
            rotcube = None
        best_ttm, max_tran = _analyze_slit_scan(datacube=rotcube, xpts=xpts, fibers=fibers, Path=Path)
        best_ttms[i] = best_ttm
        max_trans[i] = max_tran

    if rot_pts.size > 1 and plot:
        best_rot = rot_pts[np.argmax(max_trans)]
        plt.plot(rot_pts, max_trans)
        plt.axvline(x=best_rot, color='k', linestyle='--', label='Best rotator position')
        plt.xlabel('NIRSPEC Rotator Position [deg]')
        plt.ylabel('Maximum Total Slit Throughput')
        plt.title('Slit Throughput vs Rotator Position')
        plt.savefig(Path[0]+Path[1]+savename+'_rot_scan.png',bbox_inches='tight')
        plt.show()
        plt.close()

    np.savez('ttm_rotator.npz', best_ttms=best_ttms, rot_pts=rot_pts, trans=max_trans)

    ### Return best TTM positions. If no rotator specified, just a float for TTM, otherwise
    ### Tuple of best TTM, best rotator
    if best_ttms.size == 1:
        best_ttm = np.round(best_ttms[0],3)
        return best_ttm
    else:
        best_ttm = np.round(best_ttms[np.argmax(max_trans)],3)
        return best_ttm, best_rot

def slit_rot_scan(start=-0.35, stop=-0.25, step=0.01, fibers=[1,2,3,4], rstart=None, rstop=None, rstep=None, savename='', Path=None, plot=True):
    assert start < stop 
    if rstart is not None:
        assert rstart < rstop
        assert rstep is not None
    if bkgd is None:
        bkgd, _, _, = _take_backgrounds(itime, coadds, nframes=5)
    if Path is None:
        Path = Organization.get_path('Spec')
        Path = list(Path)
        Path[1] =Path[1].split('.')[0]
    if fibers is None:
        fibers = [-1]

    _,_, datacubes = _take_slit_rot_scan(start=start, stop=stop, step=step, fibers=fibers, 
                                         rstart=rstart, rstop=rstop, rstep=restep, savename=savename, Path=Path)

    best_ttm, best_rot = _analyze_slit_rot_scan(datacubes=datacubes, start=start, stop=stop, step=step, fibers=fibers, 
                                                rstart=rstart, rstop=rstop, rstep=restep, savename=savename, Path=Path, plot=plot)
    print('Best TTM positions:', best_ttm)
    print('Best rotator positions:', best_rot)
    return best_ttm, best_rot

def flux_extract(fluxbase=None, low=None, high=None, darkind=None, fullspec=False, plot=False, save=False, plotTrace=False):
    '''
    Flux extraction/integrated flux calculator for use during the night. Currently assumes a single fiber
    '''
    if fluxbase is None or low is None or high is None:
        fname = spec.get('lastfile')
        low, high = 0, 1
    else:
        fname = None

    inds = range(low, high)
    datacube = np.empty((len(inds), 2048, 2048))

    
    # dark = fits.open(fluxbase+str(darkind).zfill(4)+'.fits')[0].data.T[:,::-1]
    if darkind is None:
        dark = np.zeros((2048,2048)),
        smoothed_thermal_noise = np.zeros((2048,2048))
    else:
        bkg_names = [fluxbase+str(darkind).zfill(4)+'.fits']
        dark, smoothed_thermal_noise, badpixmap = background.make_badpixmap(bkg_names)
    badpixmap = np.ones((2048,2048))
    badpixcube = np.expand_dims(badpixmap, axis=0)

    if fname is None:
        for i, ind in enumerate(inds):
            datacube[i] = fits.open(fluxbase+str(ind).zfill(4)+'.fits')[0].data.T[:,::-1] - dark
    else:
        datacube[0] = fits.open(fname)[0].data.T[:,::-1] - dark

    ### Fast flux extraction
    ### Rectify
    rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber)
    ### Trace fitting
    trace_locs = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)
    # if trace_locs.ndim == 2:
    #   trace_locs = np.expand_dims(trace_locs,)
    ### Flux extraction
    _, flux_lst = _extract_fluxes_fast(datacube, trace_locs, plot=plotTrace)
    
    ### Plot output fluxes
    if plot:
        for i in range(len(flux_lst[0])):
            plt.plot(flux_lst[0,i], label='order'+str(i))
        plt.legend()
        plt.show()

    if save:
        print(flux_lst.shape)
        Path = Organization.get_path()
        Path = list(Path)
        if fname is not None:
            hdu = fits.PrimaryHDU(flux_lst)
            hdu.writeto(Path[0]+fname+'_fluxes.fits', overwrite=True)
        else:
            for ind in inds:
                hdu = fits.PrimaryHDU(flux_lst)
                hdu.writeto(Path[0]+fluxbase+str(ind).zfill(4)+'_fluxes.fits', overwrite=True)

    if not fullspec:
        int_flux = np.empty(flux_lst.shape[0])
        for i, frame in enumerate(flux_lst):
            int_flux[i] = np.nansum(np.nansum(frame, axis=1), axis=0)
        print('Total counts in top '+str(N_order)+' orders of spec, by frame:', np.round(int_flux,0))
        return int_flux
    else:
        # print(flux_lst.shape)
        return  flux_lst

# def calculate_throughput(data, k_mag, exptime, bb_temp=5000, plot=False):
#     """
#     Roughly estimatels throughput of data. Currently only works for K-band for one particular grating configuration!!!
#     Args:
#         data (np.array): extracted fluxes with shortest wavelengths coming first. Dimensions are (Norders x Nchannels)
#         k_mag (float): 2MASS K-band magnitude
#         exptime (seconds): exposure time of frame
#         bb_temp (float): optional, blackbody temperature assumed for stellar model
#     Returns
#         throughout (float): the 95% highest throughput calculated. Nearly the peak throughput
#     """

#     bb = models.BlackBody(temperature=bb_temp*u.K)
#     star_wv = order_wvs
#     star_model = bb(order_wvs * u.um).to(u.W/u.cm**2/u.um/u.sr, equivalencies=u.spectral_density(order_wvs * u.um)).value
#     # print(star_model)
#     guess_start = k_filt['wv'][0]
#     guess_end = k_filt['wv'][-1]
#     star_model_filt = bb(k_filt['wv'] * u.um).to(u.W/u.cm**2/u.um/u.sr, equivalencies=u.spectral_density(k_filt['wv'] * u.um)).value

#     k_zpt = 4.283E-14 # W / cm^2 / micron
#     k_flux = k_zpt * 10**(-k_mag/2.5)
#     integral = np.sum(star_model_filt * k_filt['trans'])/np.sum(k_filt['trans'])
#     norm = k_flux/integral
#     photon_energy = 6.626068e-34 * 299792458 / (star_wv * 1e-6) # Joules
#     tele_size = 76 * (100)**2 # cm^2
#     model_photonrate = star_model * norm / photon_energy * (tele_size)

#     throughputs = []
#     if plot:
#         plt.figure()
#     for wv_soln, order in zip(wv_solns, data):
#         xcoords = np.arange(order.shape[0])
#         wvs = np.poly1d(wv_soln)(xcoords)
#         dlam = wvs - np.roll(wvs, 1)
#         dlam[0] = wvs[1] - wvs[0]
        
#         model_photonrate_order = np.interp(wvs, star_wv.ravel(), model_photonrate.ravel())
#         model_photonrate_order *= exptime * dlam
#         data_photons = order * gain
#         throughput = medfilt(data_photons/model_photonrate_order,kernel_size=10)
#         throughputs.append(throughput)
#         if plot:
#             plt.plot(wvs, (throughput), 'b-')
#     if plot:
#         plt.show()
#     throughputs = np.array(throughputs)

#     return np.nanpercentile(throughputs, 95)

# def frame_throughput(fluxbase='/nfiudata/210704/Spec/Raw_Frames/nspec210704_', fnum=None, darkind=None ,kmag=5.4, exptime=60., bb_temp=5000, plot=True, plotTrace=False):
#     if fluxbase is None or fnum is None:
#         fname = spec.get('lastfile')
#         fnum = int(fname.split('_')[1][:-5])
#     low = fnum
#     if fnum is not None:
#         high = fnum+1
#     else: high = None
#     fluxes = flux_extract(fluxbase=fluxbase, low=low, high=high, darkind=darkind, fullspec=True, plotTrace=plotTrace)[0,0]
#     throughputs = calculate_throughput(fluxes, kmag, exptime, bb_temp, plot=plot)
#     return throughputs

def get_flux_ratio(plotTrace=False):
    ### Make sure we have backgrounds, then take a spec image and subtract
    bkgd, smoothed_thermal_noise, badpixmap = _take_backgrounds(itime, coadds, nframes=3)
    spec.Take_Raw_Image()
    datacube = np.expand_dim(spec.im[0].T[:,::-1] - bkgd, axis=0)

    ### Get total counts on NIRSPEC
    ### Rectify frame
    rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber)
    ### Trace fitting
    trace_locs = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)
    ### Flux extraction
    _, flux_lst = _extract_fluxes_fast(datacube, trace_locs, plot=plotTrace)
    spec_cnts = np.nansum(np.nansum(flux_lst))
    spec_flux = spec_cnts/itime
   
    ### Get Cred2 flux
    cred_flux = Acquisition.get_cred2_flux(fiber=track.get_goal()[0],nb_im=100)

    return spec_flux, cred_flux, cred_flux/spec_flux


def plot_spectra(fnames=None, plotTrace=False):
    '''
    For local testing use only, uses hardcoded filenames
    '''
    fnames = ['arcs/nspec210528_0938.fits','arcs/nspec210528_0939.fits','arcs/nspec210528_0946.fits','arcs/nspec210528_0947.fits']
    if fnames is None:
        fnames = []
        nframes = input('Enter number of raw frames >> ')
        try: nframes = int(nframes)
        except:
            print('Invalid number of frames')
            return -1
        for i in range(nframes):
            fname = input('Enter filename '+str(i+1)+'/'+str(nframes)+' >> ').replace('\'', '').replace(' ', '')
            if not os.path.exists(fname):
                print('Invalid datacube filename!')
                print(fname)
                return -1
            else:
                fnames.append(fname)
    else:
        for i, name in enumerate(fnames):
            fnames[i] = name.replace('\'', '').replace(' ', '')


    ### There needs to be something global keeping trace of this
    bkg_names = ['arcs/nspec210528_0971.fits']
    # bkg_names = ['arcs/nspec210528_0955.fits']
    dark, smoothed_thermal_noise, badpixmap = background.make_badpixmap(bkg_names)
    badpixmap = np.ones((2048,2048))
    # dark = np.zeros((2048,2048))
    badpixcube = np.expand_dims(badpixmap, axis=0)

    datacube = np.empty((len(fnames),2048,2048))
    for i, fname in enumerate(fnames):
        dataarr = fits.open(fname)[0].data
        if dataarr.ndim == 3:
            dataarr = dataarr[:,:,::-1]
            for j in range(len(dataarr)):
                dataarr[j] = dataarr[j]-dark
            datacube = dataarr
        else:
            datacube[i] = dataarr.T[:,::-1]-dark
        # goodim = datacube[i]
        # goodim[np.isnan(goodim)] = 0.
        # plt.imshow(goodim, norm=LogNorm(vmin=1, vmax=1000))
        # plt.show()

    ### DRP way
    # print('DRP')
    # trace_sigs, trace_locs, flux_lst = _reduce_drp_byframe(datacube, smoothed_thermal_noise, badpixmap,  N_order, N_fiber, stack=True)
    # flux_lst_drp = np.nansum(flux_lst, axis=0)[0]

    print('Trace Finding')
    rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber)
    ### Trace fitting
    trace_locs = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)
    if N_fiber == 1:
        trace_locs = np.expand_dims(trace_locs, axis=0)
    
    # print('Hybrid')
    # flux_lst_hybrid = np.nansum(_extract_fluxes_hybrid(datacube, smoothed_thermal_noise, badpixmap, trace_locs),axis=0)[0]
    # print(flux_lst_hybrid.shape)

    ### Flux extraction
    print('Box')
    flux_lst, cnt_lst = _extract_fluxes_fast(datacube, trace_locs, plot=plotTrace)
    print(cnt_lst.shape)
    flux_lst = np.nansum(flux_lst,axis=0)[0]
    # cnt_lst = np.nansum(cnt_lst, axis=0)[0]
    print(cnt_lst.shape)
    np.save('fiber_finding_fluxes.npy', cnt_lst)
    cnt_lst = cnt_lst[0,0,::-1,::-1]
    goodordwvs = order_wvs[::-1,::-1]

    ### Plot output fluxes - assuming only one fiber
    # for j, order in enumerate(cnt_lst[0]):
    #   plt.plot(order, label='Order '+str(j))
    #   plt.legend()
    #   plt.show()
    # for j, order in enumerate(flux_lst):
    #   plt.plot(order, label='Order '+str(j))
    #   plt.legend()
    #   plt.show()
    for j in range(N_order):
        fig = plt.figure(num = 1, figsize = (9,4))
    
        gs1 = GridSpec(1, 4)
        gs1.update(left=0.05, right=0.95, top = 0.2, bottom = 0.8, hspace=0.4, wspace=0.4)

        # Create a subplot for the image
        plt.subplot(gs1[:,:2])
        plt.plot(goodordwvs[j], cnt_lst[j]/np.nanmedian(cnt_lst[j]), label='Box')
        plt.xlabel(r'Wavelength [$\mu$m]')
        plt.subplot(gs1[:,2:])
        plt.plot(cnt_lst[j]/np.nanmedian(cnt_lst[j]), label='Box')
        plt.xlabel('Pixel')
        # plt.plot(cnt_lst[j]/np.nanmedian(cnt_lst[j]))
        # plt.plot(flux_lst_hybrid[j]/np.nanmedian(flux_lst_hybrid[j]), label='Optimal')
        # plt.plot(flux_lst_drp[j]/np.nanmedian(flux_lst_drp[j]), label='DRP')
        # plt.plot(cnt_lst[j]/np.nanmedian(cnt_lst[j])-flux_lst_hybrid[j]/np.nanmedian(flux_lst_hybrid[j]), label='Box Residual')
        # plt.plot(flux_lst_hybrid[j]/np.nanmedian(flux_lst_hybrid[j])-flux_lst_drp[j]/np.nanmedian(flux_lst_drp[j]), label='Hybrid Residual')
        plt.title('Order '+str(j))
        plt.show()



def take_fiber_TTM_scan():
    '''
    Take 4D scan over fiber x/y and FEU TTM x/y. Not user-settable for now because I'm lazy
    '''
    Path = Organization.get_path('Spec')
    Path = list(Path)
    Path[1] =Path[1].split('.')[0]
    ### Take backgrounds
    bkgd, smoothed_thermal_noise, badpixmap = _take_backgrounds(itime, coadds, nframes=3)
    badpixmap = np.ones((2048,2048))
    # badpixcube = np.expand_dims(badpixmap, axis=0)
    # bkgd = np.zeros((2048, 2048))
    # smoothed_thermal_noise = np.zeros((2048, 2048))

    ### Set up scanning points
    txpts = np.arange(-0.38, -0.335, 0.01)
    typts = np.arange(-0.02, 0.025, 0.01)
    fxpts = np.arange(3.54, 3.75, 0.05)
    fypts = np.arange(5.13, 5.34, 0.05)

    datacube = np.empty((txpts.size*typts.size*fxpts.size*fypts.size, 2048, 2048))
    cnt = 0
    tipos = FEU_TTM.get_target()
    fipos = fib.get_pose()
    for i in range(txpts.size):
        for j in range(typts.size):
            FEU_TTM.set_pos([txpts[i], typts[j]], block=True)
            for k in range(fxpts.size):
                for l in range(fypts.size):
                    fib.set_pos([fxpts[k], fypts[l], ipos[-1]])
                    spec.TakeImage(Coadds=coadds, tint=itime) #, Save=False, filename='slit_scan_step'+str(i))
                    datacube[cnt] = spec.im[0].T[:,::-1] - bkgd

                    txt = '%3.0f / 625'%float(cnt)
                    sys.stdout.write('\rFrame: ' + txt)
                    sys.stdout.flush()
    print('')
    hdu = fits.PrimaryHDU(datacube)
    hdu.writeto(Path[0] + Path[1]+'_fib_ttm_scan.fits', overwrite=True)

    FEU_TTM.set_pos(tipos, block=True)
    fib.set_pos(fipos, block=True)
