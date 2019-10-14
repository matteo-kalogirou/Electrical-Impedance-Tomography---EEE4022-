function [ y ] = spread( R, S, t1, t2)
% [ y ] = spread( R, S, thresh)
% Returns a percentage representative of the size of the reconstructed
% image with respect to the another image. Analyses element data from R and
% S to determine the number of elements above a threshold value and returns
% a percentage match between reconstructed R and simualted S.
% 
% R is the element data of one img
% S is the element data of the other img
% thresh: the activation threshhold

% if(size(R)~=size(S))
%     display('Inputs must be the same size');
%     return;
% end
% 
% if(size(R,2) > size(R,1))
%     R=R';
% end
% if(size(S,2) > size(S,1))
%     S=S';
% end
% 
% 
% actr=0;
% acts=0;
% for i = 1:size(R,1)
%     if(R(i)>=thresh)
%        actr = actr+1 ;
%     end
%     if(S(i)>=thresh)
%        acts = acts+1 ;
%     end
% end


y = size(find(R.elem_data>= t1*max(R.elem_data)),1)/...
     size(find(S.elem_data>= t2*max(S.elem_data)),1);


% y = actr/acts;


