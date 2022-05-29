% FVM 2D difusion estacionaria
% -div(nu*grad(phi)) = 0
% MNFT 2021
% Juan M. Gimenez, Santiago Márquez Damián
clear variables;
close all;

addpath 'Utils';
addpath 'meshFiles';
warning("off");

%% USER PARAMETERS %%%%%%
% Caso placa cuadrangulos
% Archivo de nodos
%file_xnod = 'xnod_placa';
% Archivo de conectividades
%file_icone = 'icone_placa_quad';
% Funcion de condiciones de borde
BCs = @placa;
% Difusividad
nu = 0.5;
%% END USER PARAMETERS %%%

% generacion de las estructuras geometricas
%make_data;
filename = 'barra1D.msh';
Mesh = gen2DMesh(filename);

% patches para condiciones de borde
bf = find(Mesh.neighbour == 0);
patches = BCs(Mesh.xnod, Mesh.faces, bf);

[AD, bD] = assemble_diffusion(Mesh, patches, nu);

phi = AD\bD;

view2d_by_ele(Mesh.xnod, Mesh.icone, phi);
axis equal;

colormap jet;
colorbar;

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
