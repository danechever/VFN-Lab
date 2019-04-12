function LSQval = VFN_An_LSQfunc(X, datPSF, datDON, modPSF, modDON)
% VFN_An_LSQfunc Calculate least squares difference between model and lab data
%   
%   - The X argument should be a vector with parameters over which minimization
%       will occur.
%
%   - Takes in PSF and Donut average radial profiles as data vectors. Runs 
%       through coupling map calculation of VFN simulation to get model coupling
%       matrix given the provided parameters (fmfd) for a PSF and Donut.
%   - Calculates average radial profiles for simulations. Scales the simulation
%       results by thpt in y and pztg in x (radial axis).
%   - Does Least Squares Calculation as:
%           sum( (datVecPSF - modelPSF)^2 + (datVecDon - modelDON)^2 )
%
%   ** NOTE: Make sure that PSF and DON rad profiles are of the same length.
%       Otherwise, one will inadvertently end up weighted more than the other
%
%   - This is meant to be called iteratively by fmincon or another optimizer
%       function trying to fit the model to the data.
%
%   NOTES FOR IMPROVEMENTS:
%   - The slope of the curves at the end do not automatically match so in the
%       current LSQ calculation. As such, there would be errors in the fit
%       beyond the fit range. This could be mitigated by adding a slope term to
%       the fit for the final point.
%
%   LSQval = VFN_An_LSQfunc(X, datPSF, datDON, modPSF, modDON)
%     Calculate LSQ value via sum of results from datVecPSF and datVecDON with
%       model PSF and model Donut respectively.
%     - 'datPSF' Structure with lab data and parameters for 'PSF' (Vortex out)
%        - datPSF.radvg:    Normalized average radial profile
%        - datPSF.rax:      Radial axis (in L/D units)
%     - 'datDON' Structure with lab data and parameters for 'Donut' (Vortex in)
%        - datDON.radvg:    Normalized average radial profile
%        - datDON.rax:      Radial axis (in L/D units)
%     - 'modPSF' Structure containing model parameters for PSF
%        - modPSF.PSFv:         E-field at fiber tip
%        - modPSF.totalPower0:  Total power on fibertip
%        - modPSF.lambdaOverD:  Samples per Lambda/D
%        - modPSF.coords:       Coordinate structure
% 	     - modPSF.lambda:       Wavelength
%     - 'modDON' Structure containing model parameters for Donut
%        - modDON.PSFv:         E-field at fiber tip
%        - modDON.totalPower0:  Total power on fibertip
%        - modDON.lambdaOverD:  Samples per Lambda/D
%        - modDON.coords:       Coordinate structure
% 	     - modDON.lambda:       Wavelength
%     - 'X'         Values over which minimization will occur:
%           X(1) = 'thpt':  Factor by which to scale the models in y (coupling).
%                       This essentially accounts for lens, fiber, etc. losses
%           X(2) = 'pztg':  V/micron value to be used for the piezo actuators. 
%                       This basically scales the models in x (radial dir)
%           X(3) = 'mfdf':  MFD/F# ratio. These are degenerate parameters in
%                       building the models
%
%     Returns:
%     - 'LSQval'    Result of LSQ calculation
%  

%% Redefine variables
%-- Define minimization variables more clearly
thpt = X(1);        % Throughput scaling factor
pztg = X(2);        % piezo-gain (radial) scaling factor
mfdf = X(3);        % MFD/F# ratio

%-- Extract data from structure
%----- PSF
radvgDat_psf = datPSF.radvg;     % Average radial profile of PSF data - vector
raxDat_psf   = datPSF.rax;       % X-axis (radial) for PSF data - vector
%----- Donut
radvgDat_don = datDON.radvg;     % Average radial profile of PSF data - vector
raxDat_don   = datDON.rax;       % X-axis (radial) for PSF data - vector


%-- Define model variables clearly
% %------------ PSF model
% E_psf = modPSF.PSFv;            % E-field at fiber tip
% Pow_psf = modPSF.totalPower0;   % Total power on fibertip
% LoD_psf = modPSF.lambdaOverD;   
% CRD_psf = modPSF.coords;
% Lam_psf = modPSF.lambda;        % Wavelength
% %------------ Donut model
% E_don = modDON.PSFv;            % E-field at fiber tip
% Pow_don = modDON.totalPower0;   % Total power on fibertip
% LoD_don = modDON.lambdaOverD;   
% CRD_don = modDON.coords;
% Lam_don = modDON.lambda;        % Wavelength

%% Generate model coupling maps (using MFD/F# provided)
%-- Define fiber diameters 
    % This is where the F#/MFD ratio comes in. From the simulation code:
    %       fiberdiam = {MFD/(1.4*lambda*F#)}*1.4
    %   That equation arises from the fact that we need to rescale the MFD based
    %   on what the ideal MFD for the given F# is. The ideal MFD is provided by
    %   1.4*lambda*F# so a different MFD is some fractional amount bigger or
    %   smaller than this
    % Some algebra leads to:
    %       fibD = MFD/(lambda*F#) = mfdf/lambda
fibD_psf = mfdf/modPSF.lambda; % units of lambda/D  
fibD_don = mfdf/modDON.lambda; % units of lambda/D 

%-- Generate fiber modes
fibmod_psf = generateSMFmode_gaussian(fibD_psf*modPSF.lambdaOverD,modPSF.coords);
fibmod_don = generateSMFmode_gaussian(fibD_don*modDON.lambdaOverD,modDON.coords);

%-- Generate coupling maps
    % Note: use 8*lambda/D instead of 3* since we want to be able to vary the
    % fiber diameters by a large amount. This, however will be time-costly. Can
    % consider decreasing. 8 was chosen as:
    %       (8/(1.4*lambda*F#))*1.4 = V ---> so use 2*V to have padding
    %       (8/(1.4*0.635*3.1))*1.4 ~ 4 ---> 2*4 = 8        where 3.1=F# for original system
coup_psf = generateCouplingMap(fibmod_psf, modPSF.PSFv, modPSF.totalPower0, 8*modPSF.lambdaOverD, modPSF.coords);
coup_don = generateCouplingMap(fibmod_don, modDON.PSFv, modDON.totalPower0, 8*modDON.lambdaOverD, modDON.coords);

%% Get radial profile of model coupling maps (and rescale)
%-- Find peak (PSF) and null (Donut)
[~, ind_psf] = VFN_An_getEta_p(coup_psf);
[~, ind_don] = VFN_An_getEta_s(coup_don,2.03);

%-- Get average radial profile
cent = [ind_psf(2), ind_psf(1)];
[radvg_psf,rvec_psf] = ...
            VFN_An_radAverage(coup_psf, cent);
cent = [ind_don(2), ind_don(1)];
[radvg_don,rvec_don] = ...
            VFN_An_radAverage(coup_don, cent);
        
%-- Rescale model coupling by thpt
radvg_psf = radvg_psf*thpt;
radvg_don = radvg_don*thpt;
        
%-- Get radial axis for models
rax_psf = rvec_psf/modPSF.lambdaOverD;
rax_don = rvec_don/modDON.lambdaOverD;

%-- Rescale model at new pzt gain of data using interpolation
% Interpolate same number of points as data vectors
if length(radvgDat_psf) ~= length(radvgDat_don)
    error('Data vectors are not of equal length so weighting will be different')
end
radvg_psf = interp1(rax_psf, radvg_psf, raxDat_psf*pztg)';
radvg_don = interp1(rax_don, radvg_don, raxDat_don*pztg)';

%% Calculate LSQ value
% Since all vectors are of the same length, can use a single summation:
LSQval = sum( (radvgDat_psf - radvg_psf).^2 + (radvgDat_don - radvg_don).^2);
end