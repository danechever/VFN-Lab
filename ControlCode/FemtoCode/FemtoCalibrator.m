
for i = 6:6

    FMTO_scale = i;
    s = VFN_FMTO_setGain(s, FMTO_scale);
    readVals = nan(100,1);
    for ii = 1:length(readVals)
        readVals(ii) = mean(startForeground(s));
    end
    fprintf('%i  |  mn: %f  |  std: %f\n',FMTO_scale, mean(readVals), std(readVals))
    if mean(readVals) > 10
        break
    end
    
end

fprintf('------------\n')

FMTO_scale = 5;
s = VFN_FMTO_setGain(s, FMTO_scale);
fprintf('%i  |  %f\n', FMTO_scale, mean(mean(startForeground(s))))

