function element=make_elem(node_pattern,num_u,num_v,num_w,inc_u,inc_v,inc_w,nnx)

%
% Synopsis: element = make_elem(node_pattern,num_u,num_v,num_w,inc_u,inc_v,inc_w,nnx)
%
% Input: node_pattern > Natural node ordering for elements
% node_pattern=[ 1 2 nnx+1 nnx+2 nny+1 nny+2 nny+nnx+1 nny+nnx+2 ]; % Node pattern 1 %2 3 4 5 6 7 8 % according to me notations

% num_u > Number of control volumes in X direction.
% num_v > Number of control volumes in Y direction.
% inc_u > Increment in cell numbers in X direction
% inc_v > Increment in cell numbers in Y direction
%
% Output: element > Element connectivity matrix.


inc=[zeros(1,size(node_pattern,2))];
e=1;
element=zeros(num_u*num_v*num_w,size(node_pattern,2));

for row=1:num_v*num_u
for col=1:num_u
element(e,:)=node_pattern+inc;
inc=inc+inc_u;
e=e+1;
end
inc = row*inc_v;
if mod(e-1,num_u*num_v)==0
node_pattern = node_pattern + nnx;
end
end