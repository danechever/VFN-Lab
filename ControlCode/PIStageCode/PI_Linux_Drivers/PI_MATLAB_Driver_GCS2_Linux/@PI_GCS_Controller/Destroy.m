function Destroy( c )
%DESTROY Controller object and unload library
if(IsConnected(c))
    c = CloseConnection(c);
end
clear c;