% Close the Labjack/Femto connection

% set gain to minimum (5)
FMTO.FMTO_scale = 5;
VFN_FMTO_setGain(FMTO);

%% Close the connection
py.labjack.ljm.close(FMTO.LJ);
