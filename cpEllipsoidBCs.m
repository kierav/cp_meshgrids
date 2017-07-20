function [cpx,cpy,cpz, sdist, bdy, bdyBCs] = cpEllipsoidBCs(x,y,z, AB, cen, ax)
%CPELLIPSOID  Closest point function for oblate/prolate ellipsoid.
%   [cpx,cpy,cpz, sdist, bdy, bdyBCs] = cpEllipsoid(x,y,z, [a b]) returns the
%   closest point and distance to (x,y,z). It also returns the endpoints
%   of the two sem-principal axes. [a b] specifies semi-principal axes lengths.  
%   If omitted it takes default values.
%
%   [cpx,cpy,cpz, sdist, bdy, bdyBCs] = cpEllipsoid(x,y,z, [a b], [xc,yc,zc]) is
%   the same but centered at (xc,yc,zc)
%
%   [cpx,cpy,cpz, sdist, bdy, bdyBCs] = cpEllipsoid(x,y,z, [a b], [], 'y') is the
%   same but revolved around the y-axis.  Can also pass 'x' or 'z'.
%
%   This version makes oblate or prolate ellipsoid of revolution
%   around the x axis.  That is, so you can only specify the major
%   axis (aligned with x axis) and one minor axis which is used for
%   the other two.
%
%   Note: returns signed distance (with negative inside).

  % defaults
  if (nargin < 4) || isempty(AB)
    AB = [1.25 0.75];
  end
  if (nargin < 5) || isempty(cen)
    cen = [0 0 0];
  end
  if (nargin < 6)
    ax = 'x';
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  P = {AB(1) AB(2)};

  switch ax
    case 'x'
      [cpx, cpy, cpz, sdist] = cpSurfOfRevolution(x, y, z, @cpEllipse, 'x', P);
      bdy1 = cpx==AB(1); bdy3 = cpx==-AB(1); 
      bdy2 = cpy==AB(2); bdy4 = cpy==-AB(2);
    case 'y'
      [cpy, cpx, cpz, sdist] = cpSurfOfRevolution(y, x, z, @cpEllipse, 'x', P);
      bdy1 = cpy==AB(1); bdy3 = cpy==-AB(1); 
      bdy2 = cpz==AB(2); bdy4 = cpz==-AB(2);
    case 'z'
      [cpz, cpy, cpx, sdist] = cpSurfOfRevolution(z, y, x, @cpEllipse, 'x', P);
      bdy1 = cpz==AB(1); bdy3 = cpz==-AB(1); 
      bdy2 = cpx==AB(2); bdy4 = cpx==-AB(2);
    otherwise
      error('invalid axis of revolution');
  end

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);

  bdy = bdy1 | bdy2 | bdy3 | bdy4;
  bdyBCs = [bdy1(:),bdy2(:),bdy3(:),bdy4(:)];