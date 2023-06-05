import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.pylab as pl
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
from Nirspec_cmds import Spec_cmds
from Star_Tracker_cmds import Tracking_cmds
from Track_Cam_cmds import TC_cmds
from FAM_cmds import FAM_cmds
import Acquisition
import ktl
import Organization
import sys
import time
import warnings
import datetime
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

# print('Current save directory:', Organization.get_path()[0])
cam   = TC_cmds()
track = Tracking_cmds()
FAM  = FAM_cmds()
spec  = Spec_cmds()
pd_sensitivity = 10
nreads = 10 #100

# connect to flip mirror keyword
serv = ktl.Service("kpic")
serv["PDPONAM"].monitor()

### These are some parameters that probably don't need to change often
N_order = 7
N_fiber = 1
nskip = 4 ### Only fit every 4 points for speed when using KPIC DRP, still give good enough trace after smoothing
itime = 1500 #30*1000
coadds = 1

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
            psfx, psfy = track.get_psf_cent()
            # get target location
            goal = track.get_psf_goal()

            # if we're close enough to goal, break
            if track.get_psf_parameters()[0] and abs(psfx - goal[1]) < .5 and abs(psfy - goal[2]) < .5:
                return 1
    else:
        return 0

def _rectify(cube, badpixmap, N_order, N_fiber, bkgd=np.zeros(0)):
    '''
    Flattens and rectifies a cube for peak finding
    '''
    ### Flatten and get rid of bad pixels
    deg = 4
    ### Badpixmap caussing issues, skipping for now
    # print('cube shape', cube.shape)
    # image = np.nansum(cube, axis=0) * badpixmap
    # Median instead of sum
    if bkgd.shape[0] == 0:
        image = np.nanmedian(cube, axis=0) * badpixmap
    else:
        image = (np.nanmedian(cube, axis=0) - bkgd) * badpixmap
    # plt.imshow(image, norm=LogNorm(vmin=1,vmax=1000))
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
        pass
        #print('Rectification failed in', len(rectfails), 'columns:', np.sort(rectfails))

    rectIm[np.isnan(rectIm)] = 0.0

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
    #if np.nanmax(int_fluxes) < 100*N_order*2048:
        #print('Warining: Low counts on spec, trace fitting may be unreliable. Check source intensity!')
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

def _extract_fluxes_fast(datacube, trace_loc, N_order=7, bpmap=np.ones((2048,2048))):
    trace_sigs = np.empty((len(datacube), N_fiber, N_order, 2048))
    flux_ests = np.empty((len(datacube), N_fiber, N_order, 2048))
    tot_cnts = np.empty((len(datacube), N_fiber, N_order, 2048))
    trace_locs = np.empty((len(datacube), N_fiber, N_order, 2048))

    # TO-DO: look into parallizing this to speed up 
    for i, image in enumerate(datacube):
        # nframes = str(len(datacube))
        # txt = '%3.0f / '%(i+1) + nframes
        # sys.stdout.write('\rExtracting frame'+txt)
        # sys.stdout.flush()
        
        # apply bad pixel map
        goodim = image * bpmap
        goodim[np.isnan(goodim)] = 0.0
        trace_amp, _, trace_sig, tot_cnt =_fit_traces_fast(goodim, trace_loc, N_order, N_fiber)
        # trace_locs[i] = trace_loc
        # flux_ests[i] = trace_sig*trace_amp*np.sqrt(2*np.pi)
        # flux_ests[flux_ests>3e5] = 0.0 ### Reject hot pixel values
        # flux_ests[flux_ests<0] = 0.0 ### Reject bad pixel values\

        trace_sigs[i] = trace_sig
        tot_cnts[i] = tot_cnt

    # print('')
    return tot_cnts

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

def _one_order_fit(traces, goodim, i, j):
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

    return amps, sigs, tots

# def _fit_traces_fast_new(goodim, traces, N_order, N_fiber):
#     rawsig = np.empty((N_fiber, N_order, 2048))
#     rawamp = np.empty((N_fiber, N_order, 2048))
#     rawtot = np.empty((N_fiber, N_order, 2048))

#     if traces.ndim == 2:
#         traces = np.expand_dims(traces,axis=0)

#     for i in range(N_fiber):
#         # parallelize N_order calc
#         pool = mp.Pool(processes=N_order)
#         res_list = []
#         for j in range(N_order):
#             res_list.append(pool.apply_async(_one_order_fit, args=(traces, goodim, i, j)))
        
#         for j in range(N_order):
#             res = res_list[j].get()
#             rawamp[i,j] = res[0]
#             rawsig[i,j] = res[1]
#             rawtot[i,j] = res[2]

#         pool.close()
#     # print(rawsig)
#     ### Fourier smoothing - drop all the high-frequency stuff
#     # smoothamp = _fourier_smooth(rawamp, N_order, N_fiber)
#     smoothsig = _fourier_smooth(rawsig, N_order, N_fiber)
#     return rawamp, rawsig, smoothsig, rawtot

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
                # res = res_list[j].get()
                # inds = res[0]
                # goodinds = res[1]
                # newvals = res[2]
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

def _take_backgrounds(itime, coadds, nframes=5, retake=False, Path=None):
    '''
    Take background frames for scan
    '''
    fam_bkgd_pos = [50, 9900]

    ### Path info
    if Path is None:
        Path = Organization.get_path('Spec')
        Path = list(Path)
        Path[1]=Path[1].split('.')[0]

    # Path is actually a list, so get the string
    real_path = Path[0]
    #print(real_path)
    
    ### Check that backgrounds already exist, if they don't force a retake
    bkg_names = []
    for i in range(nframes):
        bkg_names.append(real_path + 'bkgd_'+str(i)+'_itime'+str(itime)+'.fits')
        if not os.path.isfile(bkg_names[-1]):
            retake = True
    track_flag=False
    if retake:
        # set nreads and sampmode
        set_nspec_reads(itime)
        
        ### Stop tracking and move off-slit for backgrounds
        if track.is_tracking():
            track_flag=True
            track.stop_tracking()
            # Wait for the Tracking loop to complete its last iteration
            time.sleep(2)
        ### Get initial TTM position
        TTM_CP = FAM.get_pos()
        ### Move TTM and wait for new position
        FAM.set_pos(fam_bkgd_pos, block=True)

        ### Take backgrounds
        for i in range(nframes):
            # Take backgrounds, save for DRP
            spec.TakeImage(Coadds=coadds, tint=itime, Save=True, filename='bkgd_'+str(i)+'_itime'+str(itime), checkbadframe=False)
        FAM.set_pos([TTM_CP[0], TTM_CP[1]], block=True)

    ### Return to original TTM position and restart tracking
    bkgd, smoothed_thermal_noise, badpixmap = background.make_badpixmap(bkg_names)

    if track_flag==True:
        track.start_tracking()

    return bkgd, smoothed_thermal_noise, badpixmap

# Make trace calc separate from flux extract
def _calc_trace(datacube, badpixmap, N_order, N_fiber, bkgd=np.zeros(0), plot=False):
    
    # datacube should have be numpy array of shape (Nframes, x, y)
    # Rectify images
    rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber, bkgd=bkgd)
    trace_locs = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)

    if plot:
        collapsedim = np.nansum(datacube,axis=0)
        collapsedim[np.isnan(collapsedim)] = 0.
        plt.title('Stacked cube, with traces')
        plt.imshow(collapsedim, norm=LogNorm(vmin=1,vmax=100))
        for ordloc in trace_locs:
          plt.plot(ordloc, color='r', linestyle='--')
        plt.show()

    return trace_locs

def _total_flux_frame(this_image, badpixmap, N_order, N_fiber, bkgd=np.zeros(0), loaded_Trace=np.zeros(0), which_Order=0):
    '''
    Get the total flux on spec for a single image
    '''
    
    datacube = np.zeros((1, this_image.shape[0], this_image.shape[1]))
    
    if bkgd.shape[0] == 0:
        datacube[0] = this_image
    else:
        datacube[0] = this_image - bkgd    # remove a master bkgd
    
    # print('here0')

    # Trace fitting if not given one
    if loaded_Trace.shape[0] == 0:
        # Rectify images
        rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber)
        trace_locs = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)
    else:
        trace_locs = loaded_Trace
        # print('Using input trace locations.')
 
    # Fast flux extractiong
    if N_fiber == 1:
        trace_locs = np.expand_dims(trace_locs, axis=0)
    # print('here1')
    
    # flux_lst has shape (Nframes, Nfibers, Norder, 2048)
    flux_lst =_extract_fluxes_fast(datacube, trace_locs, N_order=N_order, bpmap=badpixmap)
    
    order_flux = np.nansum(flux_lst[0][0][which_Order])

    # print('here2')
    
    # Get total flux in each frame
    int_flux = np.empty(flux_lst.shape[0])
    for i, frame in enumerate(flux_lst):
        int_flux[i] = np.nansum(np.nansum(frame, axis=2), axis=1)[0]

    return int_flux, order_flux

def _total_flux_framepath(frame_path, badpixmap, N_order, N_fiber, bkgd=np.zeros(0), loaded_Trace=np.zeros(0), which_Order=0):
    '''
    Get the total flux on spec for each image in the cube. Moved here since it's used a bunch
    '''
    
    ### Load the frame into a datacube
    hdu = fits.open(frame_path)
    this_image = hdu[0].data.T[:,::-1]
    datacube = np.zeros((1, this_image.shape[0], this_image.shape[1]))
    
    if bkgd.shape[0] == 0:
        datacube[0] = this_image
    else:
        datacube[0] = this_image - bkgd    # remove a master bkgd
    
    ### Trace fitting if not given one
    if loaded_Trace.shape[0] == 0:
        rectIm, rectCoeffs, goodPeaks = _rectify(datacube, badpixmap, N_order, N_fiber)
        trace_locs = _find_trace(rectIm, rectCoeffs, goodPeaks, N_order, N_fiber)
    else:
        trace_locs = loaded_Trace
        # print('Using input trace locations.')
 
    # print(trace_locs.shape)
    ### Fast flux extractiong
    if N_fiber == 1:
        trace_locs = np.expand_dims(trace_locs, axis=0)
        
    # flux_lst has shape (Nframes, Nfibers, Norder, 2048)
    flux_lst =_extract_fluxes_fast(datacube, trace_locs, N_order=N_order, bpmap=badpixmap)
    
    order_flux = np.nansum(flux_lst[0][0][which_Order])
    
    ### Get total flux in each frame
    int_flux = np.empty(flux_lst.shape[0])
    for i, frame in enumerate(flux_lst):
        int_flux[i] = np.nansum(np.nansum(frame, axis=2), axis=1)[0]
    ### If averaging < 100 counts/pixel, more likely to run into trace fitting issues
    #if np.nanmax(int_flux) < 100*N_order*2048:
        #print('Warining: Low counts on spec, trace fitting may be unreliable. Check source intensity!')
      
    return int_flux, order_flux, trace_locs

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





