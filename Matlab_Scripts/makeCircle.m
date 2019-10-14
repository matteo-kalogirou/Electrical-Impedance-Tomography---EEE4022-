function [ C ] = makeCircle( )
% Function generates a circle described by circle center (x,y) and radius
% r. These parameters are randomly determined and returned in a circle
% object C.

C.r = randi([10, 40], 1)/100;          % radius [0.1, 0.4]

fits = false;

while(~fits)
   C.x = randi([-100, 100], 1)/100;    % x [-1, 1]
   C.y = randi([-100, 100], 1)/100;    % x [-1, 1]
   if( sqrt(C.x^2 + C.y^2) <= (0.9-C.r) )
      fits =true; 
   end
end



end

