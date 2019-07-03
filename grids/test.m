numx = 10;
numy = numx;
x=[0:1/(numx):1];

y=[0:1/(numy):1]; %Matlab's meshgrid is used to create 2D grid from specified divisons above
[X,Y] = meshgrid(x,y);


X1=reshape(X',length(x)*length(y),1);

Y1=reshape(Y',length(x)*length(y),1); %Coordinates of the node

node=[X1 Y1]; % Node
tri = triangulate(X1,Y1,element); % element connectivity

nele = size(tri,1);

Z= zeros(length(y)-2,length(x)-2);

nn = tri; % nodes in a long list

xx = X1(nn); yy = Y1(nn); % coordinates corresponding to nodes

xplot = reshape(xx,size(tri));

yplot = reshape(yy,size(tri));figure(1);

clf;

fill(xplot',yplot','w');

title(['Triangular Mesh ', num2str(length(nn)),' elements'])

%%%%%%Subfunction % triangulate

