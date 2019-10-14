function [ isIn ] = square( x,y, ~ )
% Function determines wheher a point lies within the square that is deined
% outside of this function. A square object must have been created prior to
% this function being called.
% Function returns a boolean value;

isIn = false;

if(y <= s.l1(1)*x + s.l1(2) && y >= s.l2(1)*x + s.l2(2))
    if(y <= s.l4(1)*x + s.l4(2) && y >= s.l3(1)*x + s.l3(2))
        isIn = true;
    end
end


end

