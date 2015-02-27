function [ x ] = DivDiff( T, Y )

if ( size(T) == 1)
    x = Y;
else
    x = (DivDiff(T(2:end), Y(2:end)) - DivDiff(T(1:end-1),Y(1:end-1)))/(T(end) - T(1));
end
end

