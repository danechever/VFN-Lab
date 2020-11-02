function hasOpenDaisyChain = HasOpenDaisyChainConnection (c)

if (c.DC_ID > -1)
    hasOpenDaisyChain = true;
else
    hasOpenDaisyChain = false;
end