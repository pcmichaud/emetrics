

function tri = triangulate(xnod,ynod,nodes)
%
% tri = triangulate(xnod,ynod,nodes)
%
% triangulate the quadrilateral mesh
%

nele = size(nodes,1);
tri = zeros(3,2*nele)';
iv = [];
ii1 = [2 3 1];
jj1 = [4 1 3];
ii2 = [1 2 4];
jj2 = [2 3 4];
nrtri = 0;
for iel = 1:nele
iv = nodes(iel,:);
d1 = norm([xnod(iv(1))-xnod(iv(3));ynod(iv(1))-ynod(iv(3))]);
d2 = norm([xnod(iv(2))-xnod(iv(4));ynod(iv(2))-ynod(iv(4))]);
if d1 <= d2
nrtri = nrtri+1;
tri(nrtri,:) = iv(ii1);
nrtri = nrtri+1;
tri(nrtri,:) = iv(jj1);
else
nrtri = nrtri+1;
tri(nrtri,:) = iv(ii2);
nrtri = nrtri+1;
tri(nrtri,:) = iv(jj2);
end
end