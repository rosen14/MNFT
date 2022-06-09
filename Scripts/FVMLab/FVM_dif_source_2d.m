% FVM 2D
% -div(nu*grad(phi)) = S
% MNFT 2022
% Rosendo Jesús Fazzari
clear variables;
close all;

addpath 'Utils';
addpath 'meshFiles';
warning("off");

%% USER PARAMETERS %%%%%%
% Funcion de condiciones de borde
BCs = @placa;
% Difusividad
nu = 0;
Q = 0;
deltaT = 0.01;
t_max = 0.5;
theta = 0.5;

%% END USER PARAMETERS %%%

% generacion de las estructuras geometricas
%make_data;
filename = 'barra1D.msh';
Mesh = gen2DMesh(filename);

% patches para condiciones de borde
bf = find(Mesh.neighbour == 0);
patches = BCs(Mesh.xnod, Mesh.faces, bf);

[Ad, bd] = assemble_diffusion(Mesh, patches, nu);
bs = assemble_source(Mesh, patches, Q);

% Construyo el array de velocidad (un vector para cada cara)
%v = zeros(Mesh.nfaces, 2);
%v(:, 1) = 1;
v = ones(Mesh.nfaces, 2);
% Término advectivo
[Aa, ba] = assemble_advection(Mesh, patches, v, 'CD');

A = Ad + Aa;
b = bd + bs + ba;


M = Mesh.V.*eye(Mesh.ncells)/deltaT + A*theta;
M_inv = inv(M);

t = 0;
phi_init = zeros(Mesh.ncells,1);
phi_k = phi_init;
while t < t_max
  R = b + phi_k.*Mesh.V/deltaT - A*(1 - theta)*phi_k;
  phi_k1 = M_inv*R;
  t = t + deltaT;
  phi_k = phi_k1;
end

phi = A\b;

view2d_by_ele(Mesh.xnod, Mesh.icone, phi_k);
axis equal;

colormap jet;
colorbar;

figure(2);clf;plot(Mesh.xnod(24:42),phi_k(2:20))
#phi(1:10)

%  %ejemplo postproceso sampling sobre curvas en mallas de triangulos
%  TR = triangulation(icone,xnod);
%  yy=0:0.05:1;xx=zeros(size(yy));
%  ic = TR.pointLocation(xx',yy');
%  figure(2);clf;plot(yy,phi(ic));legend('phi sobre el eje central')

% % ejemplo postproceso general
% yy=0:0.05:1;xx=zeros(size(yy));
% phin = ele_field_to_nod_field(xnod,icone,phi);
% phiC = griddata(xnod(:,1),xnod(:,2),phin,xx,yy);
% figure(2);clf;plot(yy,phiC);legend('phi sobre el eje central')


%Mesh.xnod(24:42)
%figure(2);clf;plot(Mesh.xnod(24:42),phi(2:19))
