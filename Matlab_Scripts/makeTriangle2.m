function [ T ] = makeTriangle2( name )
% Creates a triangle object by randomly selecting three points, rotating
% them and then placing them within the bounds of a unit circle. The
% equaiton describing the triangle is stored in TraingleObject.fcn
%   Supply the name of the object to correctly name the output string.

% Cite the InterX function used

%  Default name  = tri
if isempty(name)
    name = 'tri';
end
    
T.l = randi([10, 50], 1)/100;
T.theta = randi([1, 90],1)/100;

% Rotation transformation matrix
T.R = [cosd(T.theta) -sind(T.theta); sind(T.theta) cosd(T.theta)];

space = 0.2;
fits = false;

while(~fits)

    T.p(1).x = randi([fix(-T.l*100), fix(T.l*100)], 1)/100;
    T.p(1).y = randi([fix(-T.l*100), fix(T.l*100)], 1)/100;

    T.p(2).x = randi([fix(-T.l*100), fix(T.l*100)], 1)/100;
    T.p(2).y = randi([fix(-T.l*100), fix(T.l*100)], 1)/100;

    T.p(3).x = randi([fix(-T.l*100), fix(T.l*100)], 1)/100;
    T.p(3).y = randi([fix(-T.l*100), fix(T.l*100)], 1)/100;
    
    if( dist(T.p(1).x, T.p(1).y, T.p(2).x, T.p(2).y) >=space && ...
        dist(T.p(1).x, T.p(1).y, T.p(3).x, T.p(3).y) >=space && ...
        dist(T.p(2).x, T.p(2).y, T.p(3).x, T.p(3).y) >= space )
        fits =true;
    end
  
end



% Apply rotation
for i = 1:3
    r = T.R*[T.p(i).x; T.p(i).y];     % Apply rotation
    T.p(i).rot_x = r(1); T.p(i).rot_y = r(2);
end

% Apply Offset
max_shift = 0.9- T.l;
fits = false;
while(~fits)
   x_shift = randi([ fix(-100*max_shift), fix(100*max_shift)] )/100;
   y_shift = randi([ fix(-100*max_shift), fix(100*max_shift)] )/100;

    for i = 1:3 
        T.p(i).trans_x = T.p(i).rot_x + x_shift;
        T.p(i).trans_y = T.p(i).rot_y + y_shift;
    end
    
    if( dist(T.p(1).trans_x, T.p(1).trans_y, 0, 0) < 0.95 && ...
        dist(T.p(1).trans_x, T.p(1).trans_y, 0, 0) < 0.95 && ...
        dist(T.p(2).trans_x, T.p(2).trans_y, 0, 0) < 0.95        )
        fits = true;
    end

end


% Characterise triangle
T.l1 = polyfit([T.p(2).trans_x, T.p(1).trans_x], [T.p(2).trans_y, T.p(1).trans_y], 1);
T.l2 = polyfit([T.p(1).trans_x, T.p(3).trans_x], [T.p(1).trans_y, T.p(3).trans_y], 1);
T.l3 = polyfit([T.p(2).trans_x, T.p(3).trans_x], [T.p(2).trans_y, T.p(3).trans_y], 1);

if( T.p(3).trans_y >= T.l1(1)*T.p(3).trans_x + T.l1(2) )
    T.eq1 =  '>=';
else
    T.eq1 =  '<=';
end
if( T.p(2).trans_y >= T.l2(1)*T.p(2).trans_x + T.l2(2) )
    T.eq2 =  '>=';
else
    T.eq2 =  '<=';
end
if( T.p(1).trans_y >= T.l3(1)*T.p(1).trans_x + T.l3(2) )
    T.eq3 =  '>=';
else
    T.eq3 =  '<=';
end

T.fcn = [ '@(x,y,z) ( y ' T.eq1 'x*' name '.l1(1) +' name '.l1(2) & ' ...
                   'y ' T.eq2 'x*' name '.l2(1) +' name '.l2(2) & ' ...
                   'y ' T.eq3 'x*' name '.l3(1) +' name '.l3(2) )' ];
               
% % Find line intercepts

x1 = linspace(-1,1,101);
y1 = T.l1(1)*x1 + T.l1(2);
y2 = T.l2(1)*x1 + T.l2(2);
y3 = T.l3(1)*x1 + T.l3(2);

%Intersection of l1 with l2 and l3
T.intX_l1 = [InterX([x1;y1], [x1;y2]), InterX([x1;y1], [x1;y3])];
T.intX_l2 = [InterX([x1;y2], [x1;y1]), InterX([x1;y2], [x1;y3])];
T.intX_l3 = [InterX([x1;y3], [x1;y1]), InterX([x1;y3], [x1;y2])];



end

