function [ output ] = zeroCornerMatrix( input )
    % Makes the four corners = 0
    output = input;
    output(1,1)     = 0;
    output(1,end)   = 0;
    output(end,1)   = 0;
    output(end,end) = 0;
end
