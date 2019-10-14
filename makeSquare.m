function [ S ] = makeSquare()
% Creates a square object, S, with randomized side length, angle of
% rotation and offset from the centre.
% Each of the sides is described by a 1st order polynomial in the form:
%           l = mx+c
% Which are described by the supporting vertices, p(1-4).
% The square is constructed by first randomizing a side length and placing
% the vertices centred around the origin.
% Next the poitns are rotated through [0;90] degrees and a randomized
% offset is applied.
% Finally the polynomials describing the lines are defined and the square
% object is returned.

S.length = randi([10, 50], 1)/100;           % length [0.1, 0.5]
S.theta = randi([1,89],1);                   % theta [1, 89] 

% Rotation transformation matrix
S.R = [cosd(S.theta) -sind(S.theta); sind(S.theta) cosd(S.theta)];

% Place points
S.p(1).x = 0-S.length/2;             %top left
S.p(1).y = 0+S.length/2;

S.p(2).x = 0+S.length/2;             %top right
S.p(2).y = 0+S.length/2;

S.p(3).x = 0+S.length/2;             %bottom right
S.p(3).y = 0-S.length/2;

S.p(4).x = 0-S.length/2;             %bottom left
S.p(4).y = 0-S.length/2;

% Apply Rotation
for i = 1:4
    r = S.R*[S.p(i).x; S.p(i).y];     % Apply rotation
    S.p(i).rot_x = r(1); S.p(i).rot_y = r(2);
end

% Apply Offset
max_shift_dist = 0.95-sqrt(S.p(1).x^2 + S.p(1).y^2);
fits = false;
while(~fits)
   x_shift = randi([-100, 100])/100;
   y_shift = randi([-100, 100])/100;
   if(sqrt(x_shift^2 + y_shift^2) <= max_shift_dist)
       fits = true;
   end
end

for i = 1:4  
    S.p(i).trans_x = S.p(i).rot_x + x_shift;
    S.p(i).trans_y = S.p(i).rot_y + y_shift;
end

% Characterise square
% l1 = p1 -> p2 (l1||l2) l2 = p4 -> p3
% l3 = p1 -> p4 (l3||l4) l4 = p2 -> p3
S.l1 = polyfit([S.p(1).trans_x, S.p(2).trans_x], [S.p(1).trans_y, S.p(2).trans_y], 1);
S.l2 = polyfit([S.p(4).trans_x, S.p(3).trans_x], [S.p(4).trans_y, S.p(3).trans_y], 1);
S.l3 = polyfit([S.p(1).trans_x, S.p(4).trans_x], [S.p(1).trans_y, S.p(4).trans_y], 1);
S.l4 = polyfit([S.p(2).trans_x, S.p(3).trans_x], [S.p(2).trans_y, S.p(3).trans_y], 1);

% Find itnercept points
x1 = linspace(-1,1,101);
y1 = S.l1(1)*x1 + S.l1(2);
y2 = S.l2(1)*x1 + S.l2(2);
y3 = S.l3(1)*x1 + S.l3(2);
y4 = S.l4(1)*x1 + S.l4(2);

%Intersection of l1 with l2 and l3
S.intX_l1 = [InterX([x1;y1], [x1;y3]), InterX([x1;y1], [x1;y4])];
S.intX_l2 = [InterX([x1;y2], [x1;y3]), InterX([x1;y2], [x1;y4])];
S.intX_l3 = [InterX([x1;y3], [x1;y1]), InterX([x1;y3], [x1;y2])];
S.intX_l4 = [InterX([x1;y4], [x1;y1]), InterX([x1;y4], [x1;y2])];

end

