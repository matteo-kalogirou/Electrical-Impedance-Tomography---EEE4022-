function [ pltx, plty ] = makeCirclePlot( x,y,r )
% Returns the x and y vectors used to plot a circle
% x = x circle centre
% y = y circle centre
% r = radius

theta = 0 : 0.01 : 2*pi;
pltx = r * cos(theta) + x;
plty = r * sin(theta) + y;

end
