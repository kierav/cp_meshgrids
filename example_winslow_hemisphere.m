%% Winslow variable diffusion mesh generator problem on hemisphere
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example generates mesh grids on a hemisphere according to Winslow
% variable diffusion.


dx = 0.2/2^1;  % grid spacing
R=1; %Radius of the circle

cpf = @cpHemisphereBCs;  
paramf = @paramHemisphere;
loadData = 1;

pad = 5;
x1d=(-R-pad*dx):dx:(R+pad*dx);
y1d=(-R-pad*dx):dx:(R+pad*dx);
z1d=(-R-pad*dx):dx:(R+pad*dx);
nx=length(x1d);
ny=length(y1d);
nz=length(z1d);
[x3d,y3d,z3d]=meshgrid(x1d,y1d,z1d);

% using the standard CP function, we get a homogeneous Neuamann BC
% [Ruuth & Merriman 2008]
%[cpx,cpy,cpz, dist, bdy] = cpHemisphere(x3d,y3d,z3d R);

% Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:
% [cpx,cpy,cpz, dist, bdy, bdyBCs] = cpf(x3d,y3d,z3d, R);
[cpxb,cpyb,cpzb,dist,bdy,bdyBCs] = cpbar_3dBCs(x3d,y3d,z3d,cpf,R);

dim = 3;
p = 3;  % degree interp
order = 2;  % laplacian order
bw = 1.001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
if (pad < bw)
  warning('maybe insufficient padding');
end
band = find(abs(dist) <= bw*dx);
outband = find(abs(dist) > bw*dx);

% Last place we need the 3D array
xg = x3d(band);
yg = y3d(band);
zg = z3d(band);
% cpxg = cpx(band);
% cpyg = cpy(band);
% cpzg = cpz(band);
bdyg = bdy(band);
cpxgb = cpxb(band);
cpygb = cpyb(band);
cpzgb = cpzb(band);
bdy1 = bdyBCs(:,1);
bdy1 = bdy1(band);
bdy2 = bdyBCs(:,2);
bdy2 = bdy2(band);
bdy3 = bdyBCs(:,3);
bdy3 = bdy3(band);
bdy4 = bdyBCs(:,4);
bdy4 = bdy4(band);


%% discrete operators
disp('building laplacian and interp matrices');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band,band);
Eu = interp3_matrix(x1d,y1d,z1d,cpxgb,cpygb,cpzgb, p, band);
Ev = Eu;
I = speye(size(Eu));

% Dirichlet BCs: mirror for ghost points outside of surface edges.
% Comment this out for Neumann BCs.
Eu(bdy1,:) = -Eu(bdy1,:);
Eu(bdy3,:) = -Eu(bdy3,:);
Ev(bdy2,:) = -Ev(bdy2,:);
Ev(bdy4,:) = -Ev(bdy4,:);

%% plotting grid
[xp,yp,zp] = paramf(48, R);
xp1 = xp(:);  yp1 = yp(:);  zp1 = zp(:);
disp('building plotting matrices');
Eplot = interp3_matrix(x1d,y1d,z1d, xp1,yp1,zp1, p, band);

figure(1); clf;


%% Do some calculation
% timestepping or whatever using the cp_matrices.

%% parameters and functions for Winslow
%   $$u_t = Del * (rho_u(x,y)*Del(u))$$
%   $$v_t = Del * (rho_v(x,y)*Del(v))$$
rho_u = @(x,y) 1;
rho_v = @(x,y) 1;


%% initial conditions - perturbation from steady state
[thg,rg,zg] = cart2pol(xg,yg,zg);
pert = 0.5*exp(-(10*(zg-.1)).^2) + 0.5*rand(size(xg));
u0 = 1 - pert;  v0 = 0.5*pert;
u0(thg >= 0 & thg <= pi/2) = 0;
u0(thg >= -pi & thg <= -pi/2) = 1;
v0(thg >= pi/2 & thg <= pi) = 0;
v0(thg >= -pi/2 & thg <= 0) = 1;
u = u0;  v = v0;

%% initial boundary conditions
u(bdy1) = 0;
u(bdy3) = 1;
v(bdy2) = 0;
v(bdy4) = 1;

%% time-stepping
Tf = 3;
dt = 1/6 * dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;

figure(1);
sphplot = Eplot*u;
sphplot = reshape(sphplot, size(xp));
Hplot = surf(xp, yp, zp, sphplot);
title('initial u')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
view(-10, 60)
shading interp
camlight left
colorbar
pause

%% Method-of-lines approach
% See [vonGlehn/Macdonald/Maerz 2013]
%lambda = 6*max(nuu,nuv)/(dx^2);
%Au = nuu*(E*L) - lambda*(I-E);
%Av = nuv*(E*L) - lambda*(I-E);

for kt = 1:numtimesteps
  %% MOL: explicit Euler timestepping
  %unew = u + dt*( E*f(u,v) + Au*u );
  %vnew = v + dt*( E*g(u,v) + Av*v );
  %u = unew;
  %v = vnew;
  %% MOL: without precomputing matrices
  %rhsu = nuu*(L*u) + f(u,v);
  %rhsv = nuv*(L*v) + g(u,v);
  %unew = u + dt*( E*rhsu - lambda*(u - E*u) );
  %vnew = v + dt*( E*rhsv - lambda*(v - E*v) );
  %u = unew;
  %v = vnew;

  %% Ruuth-Merriman
  unew = u + dt*(L*u);
  vnew = v + dt*(L*v);
  u = Eu*unew;
  u(bdy3) = 2 + u(bdy3);
  v = Ev*vnew;
  v(bdy4) = 2 + v(bdy4);

  t = kt*dt;

  if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
    disp([kt t]);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    drawnow;
  end
end

%% plot meshgrid
Nu = 256;
Nv = 256;
% Egrid = interp3_matrix(x1d,y1d,z1d,xg,yg,zg,p);
% ugrid = Egrid\u;
% ugrid = reshape(ugrid,size(x3d));
% vgrid = Egrid\v;
% vgrid = reshape(vgrid,size(x3d));
ugrid = -2*ones(size(x3d));
ugrid(band) = u;
vgrid = -2*ones(size(x3d));
vgrid(band) = v;

figure
hold on
contourslice(x1d,y1d,z1d,ugrid,xp,yp,zp,Nu)
contourslice(x1d,y1d,z1d,vgrid,xp,yp,zp,Nv)
view(3)

%% determine grid points
tol = 0.002;
[xp_,yp_,zp_] = paramf(1024,R);
xp_ = xp_(:); yp_ = yp_(:); zp_ = zp_(:);
Emesh = interp3_matrix(x1d,y1d,z1d,xp_,yp_,zp_,p,band);
uplt = Emesh*u;
vplt = Emesh*v;

M = zeros(Nu*Nv,5);
% Iu = sparse(size(uplt),Nu);
% Iv = sparse(size(vplt),Nv);
for i = 1:Nu
    Iu = abs(uplt-(i)/(Nu+1))<tol;
    for j = 1:Nv
        Iv = abs(vplt-(j)/(Nv+1))<tol;
        I = Iu&Iv;
        if sum(I) ~= 0
            M((i-1)*Nv+j,:) = [mean(xp_(I)),mean(yp_(I)),mean(zp_(I)),mean(uplt(I)),mean(vplt(I))];
        else
            M((i-1)*Nv+j,:) = [0,0,0,(i)/(Nu+1),(j)/(Nv+1)];
        end
    end
end

% plot3(M(:,1),M(:,2),M(:,3),'*');

%% Plot image
img = imread('m345-jkv.jpg');
imgr = reshape(img,[256*256,3]);

% figure
% scatter3(M(:,1),M(:,2),M(:,3),10,imgr/255)

% convert to grayscale
gr = rgb2gray(img);
% interpolate
uu = linspace(0,1,256);
vv = uu;
[UU,VV] = meshgrid(uu,vv);
wp = interp2(UU,VV,double(gr),uplt,vplt);

figure
surf(xp_,yp_,zp_,reshape(wp,size(xp_)));
colormap gray
axis equal; shading interp; colorbar;