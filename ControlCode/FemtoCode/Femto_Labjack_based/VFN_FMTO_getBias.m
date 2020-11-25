function bias = VFN_FMTO_getBias(FMTO)
% VFN_FMTO_GETBIAS Returns the known Femto bias at a given gain setting:
%
%   * Uses the bias values measured previously.
%
%   - FMTO.FMTO_scale must be within the bounds of the power meter:
%       [5 <= FMTO.FMTO_scale <= 11]
%     - If not, a warning is thrown and NaN is returned
%
%   - Assumes the gain setting is set to "high-speed" on the Femto
%
%   EXAMPLE:______
%   bias = VFN_FMTO_GETBIAS(FMTO_scale) 
%       FMTO:       struct of Femto parameters. Needs field:
%                   - FMTO_scale: Current scale setting defined by 10^FMTO_scale
%       bias:       Bias at this gain setting
%
%   See also: VFN_setUpFMTO, VFN_cleanUpFMTO

switch FMTO.FMTO_scale
    case 5
        bias = 0.004905;%0.0033887;
    case 6
        bias = 0.005382;%0.0042960;
    case 7
        bias = 0.005405;%0.0044690;
    case 8
        bias = 0.005350;%0.0047430;
    case 9
        bias = 0.004941;%0.0068927;
    case 10
        bias = -0.001120;
    case 11
        bias = -0.059031;
    otherwise
        warning('No bias for FMTO_scale = %i\n',FMTO.FMTO_scale)
        bias = nan;
end

end