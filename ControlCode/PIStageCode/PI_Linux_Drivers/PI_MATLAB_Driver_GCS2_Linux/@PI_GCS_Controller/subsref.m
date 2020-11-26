function varargout = subsref(c, S)
% function subsref(c, S)
% MATLAB uses the built-in subsref function to interpret indexed references to objects. To modify the indexed reference behavior of objects, overload subsref in the class.
% see: http://www.mathworks.com/help/matlab/ref/subsref.html

% This code is provided by Physik Instrumente(PI) GmbH&Co.KG
% You may alter it corresponding to your needs
% Comments and Corrections are very welcome
% Please contact us by mailing to support-software@pi.ws

par = '';
for n = 1:length(S)
    switch S(n).type
        case '.'
            fun = S(n).subs;
        case '()'
            par = S(n).subs;
    end
end


outsize = nargout;

% Handles the number of outputs
% "helpMessage = PIdevice.qHLP()" will alway return the message (because narout == 1), but 
% "PIdevice.qHLP()" would return nothing (because nargout == 0). But the
% following code will handle this for certain commands.

if((outsize == 0)&&(~isempty([...
        regexp(fun,'q[A-Z]..'),...
        regexp(fun,'Is[A-Z]..*'),...
        regexp(fun,'Get[A-Z]..*'),...
        regexp(fun,'Enumerate..*'),...
        regexp(fun,'ListConnectedDaisyChainDevices'),...
        regexp(fun,'HasOpenDaisyChainConnection'),...
        regexp(fun,'Receive[A-Z]..*')])))
    outsize = 1;
end

if(outsize == 1)
    
    switch(length(par))
        case 0
            varargout{1} = feval(fun,c);
        case 1
            varargout{1} = feval(fun,c,par{1});
        case 2
            varargout{1} = feval(fun,c,par{1},par{2});
        case 3
            varargout{1} = feval(fun,c,par{1},par{2},par{3});
        case 4
            varargout{1} = feval(fun,c,par{1},par{2},par{3},par{4});
        case 5
            varargout{1} = feval(fun,c,par{1},par{2},par{3},par{4},par{5});
        case 6
            varargout{1} = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6});
        case 7
            varargout{1} = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7});
        case 8
            varargout{1} = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8});
        case 9
            varargout{1} = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8},par{9});
        case 10
            varargout{1} = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8},par{9},par{10});
    end
    
elseif(outsize == 2)
    switch(length(par))
        case 0
            [varargout{1},varargout{2}] = feval(fun,c);
        case 1
            [varargout{1},varargout{2}] = feval(fun,c,par{1});
        case 2
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2});
        case 3
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3});
        case 4
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3},par{4});
        case 5
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3},par{4},par{5});
        case 6
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6});
        case 7
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7});
        case 8
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8});
        case 9
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8},par{9});
        case 10
            [varargout{1},varargout{2}] = feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8},par{9},par{10});
    end
    
elseif(outsize == 0)
    switch(length(par))
        case 0
            feval(fun,c);
        case 1
            feval(fun,c,par{1});
        case 2
            feval(fun,c,par{1},par{2});
        case 3
            feval(fun,c,par{1},par{2},par{3});
        case 4
            feval(fun,c,par{1},par{2},par{3},par{4});
        case 5
            feval(fun,c,par{1},par{2},par{3},par{4},par{5});
        case 6
            feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6});
        case 7
            feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7});
        case 8
            feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8});
        case 9
            feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8},par{9});
        case 10
            feval(fun,c,par{1},par{2},par{3},par{4},par{5},par{6},par{7},par{8},par{9},par{10});
    end
end