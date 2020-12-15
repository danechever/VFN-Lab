function [V] = VFN_FMTO_readN(FMTO)
% VFN function to read the Femto output a given number of times
%
%    Args: FMTO = a struct representing the Femto object and its params
%    Returns: V = a vector of length FMTO.Nread containing read values

    V = zeros(1,FMTO.Nread);
    for II=1:FMTO.Nread
        V(II) = py.labjack.ljm.eReadName(FMTO.LJ, FMTO.AINPort);
    end
end