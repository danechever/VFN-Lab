function VFN_NKT_setWvlRange(NKT,lower_nm,upper_nm)
%   Inputs:
%       'NKT' - object containing NKT-related stuff
%       'lower_nm' - lower wavelength (nm)
%       'upper_nm' - upper wavelength (nm)
%
%   author: G. Ruane
%   last modified: March 26, 2019
%   last modified: D. Echeverri June 11, 2019

    NKT_MIN_WVL = 400;% nm
    NKT_MAX_WVL = 900;% nm

    % Get the current lower and upper wavelengths  
    [currentlower,currentupper] = VFN_NKT_getWvlRange(NKT);
    
	% Check that the wavelengths make sense, if not, keep them the same
    if(or(lower_nm < NKT_MIN_WVL,lower_nm > NKT_MAX_WVL))
        lower_nm = currentlower;
        error('     Varia: lower wavelength outside of the allowable range.');
    end
    if(or(upper_nm < NKT_MIN_WVL,upper_nm > NKT_MAX_WVL))
        upper_nm = currentupper;
        error('     Varia: upper wavelength outside of the allowable range.');
    end

    % If the passed wavelengths are different than the current wavelengths,
    % update them
    didsomething = false;
    if(currentlower ~= lower_nm)
        NKT.nktobj.set_varia_lwpsetpoint(lower_nm);
        didsomething = true;
    end
	if(currentupper ~= upper_nm)
        NKT.nktobj.set_varia_swpsetpoint(upper_nm);
        didsomething = true;
    end

    if(didsomething)
        pause(NKT.delay);
    
%         [currentlower,currentupper] = tb_NKT_getWvlRange(bench);
%         
%         disp(['     Varia: wavelength range set to ',...
%             num2str(currentlower),'-',num2str(currentupper),'nm.']);
    end
end
