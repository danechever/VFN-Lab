function conversionGain = VFN_getFmtoResponsivity(wvl)
    % VFN_getFmtoResposivity returns the conversion gain for the femto at a
    % wavelength between 400 and 900
    % NOTE: This function is an approximation from the provided graphic
    % from femto using linear approximation with the nearest 2 points to
    % the wavelength desired
    
    % Responsivity between 400 and 900nm spaced evenly
    points = [.2258 .4139 .6188 .8049 .9282 1.0411];
    
    % Index with wavelength below wvl argument
    botInd = floor(wvl/100)-3;
    
    % Using a linear fit with the 2 nearest wavelengths recorded
    conversionGain = points(botInd)-(wvl-(botInd+3)*100)*((points(botInd)-points(botInd+1))/(100));
end