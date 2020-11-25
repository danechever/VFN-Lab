% Set gain back to 5
VFN_FMTO_LUCI_setGain(5);

% Unload LUCI library
if (libisloaded('LUCI_10_x64'))
    unloadlibrary LUCI_10_x64
end