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
deltaT = 0.001;
t_max = 0.1;
theta = 0.5;




%% END USER PARAMETERS %%%

% generacion de las estructuras geometricas
%make_data;
%filename = 'plano_triang.msh';
%filename = 'plano_triang_regular.msh';
filename = 'barra1D.msh';
Mesh = gen2DMesh(filename);

% Construyo el array de velocidad (un vector para cada cara)
v = zeros(Mesh.nfaces, 2);
v(:, 1) = 1;



% patches para condiciones de borde
bf = find(Mesh.neighbour == 0);
patches = BCs(Mesh.xnod, Mesh.faces, bf);

%phi_init = zeros(Mesh.ncells,1);

choice = 'suave';
phi_init = initial_condition(Mesh, choice);


[Ad, bd] = assemble_diffusion(Mesh, patches, nu);
bs = assemble_source(Mesh, patches, Q);


%{'upwind'; 'superbee'; 'minmod'}

tvd_sch = 'upwind';
phi_k = phi_init;
t = 0;
while t <= t_max
  t
  Mesh = interp_phif(Mesh, patches, phi_k, v, tvd_sch);
  [Aa, ba] = assemble_advection_hr(Mesh, patches, v, Mesh.psif);
  A = Ad + Aa;
  b = bd + bs + ba;

  M = Mesh.V.*eye(Mesh.ncells)/deltaT + A*theta;
  M_inv = inv(M);
  R = b + phi_k.*Mesh.V/deltaT - A*(1 - theta)*phi_k;
  phi_k1 = M_inv*R;
  t = t + deltaT;
  phi_k = phi_k1;
end
phi_k_upwind = phi_k1;

tvd_sch = 'superbee';
phi_k = phi_init;
t = 0;
while t <= t_max
  t
  Mesh = interp_phif(Mesh, patches, phi_k, v, tvd_sch);
  [Aa, ba] = assemble_advection_hr(Mesh, patches, v, Mesh.psif);
  A = Ad + Aa;
  b = bd + bs + ba;

  M = Mesh.V.*eye(Mesh.ncells)/deltaT + A*theta;
  M_inv = inv(M);
  R = b + phi_k.*Mesh.V/deltaT - A*(1 - theta)*phi_k;
  phi_k1 = M_inv*R;
  t = t + deltaT;
  phi_k = phi_k1;
end
phi_k_superbee = phi_k1;

tvd_sch = 'minmod';
phi_k = phi_init;
t = 0;
while t <= t_max
  t
  Mesh = interp_phif(Mesh, patches, phi_k, v, tvd_sch);
  [Aa, ba] = assemble_advection_hr(Mesh, patches, v, Mesh.psif);
  A = Ad + Aa;
  b = bd + bs + ba;

  M = Mesh.V.*eye(Mesh.ncells)/deltaT + A*theta;
  M_inv = inv(M);
  R = b + phi_k.*Mesh.V/deltaT - A*(1 - theta)*phi_k;
  phi_k1 = M_inv*R;
  t = t + deltaT;
  phi_k = phi_k1;
end
phi_k_minmod = phi_k1;


phi_func = zeros(Mesh.ncells,1);
if strcmp(choice, 'step')
  for cell = 1:Mesh.ncells
    if ((1/6 + t_max < Mesh.C(cell, 1)) && (Mesh.C(cell, 1) < 0.5 + t_max))
      phi_func(cell) = 1;
    endif
  endfor
else
   for cell = 1:Mesh.ncells
     if ((1/6 + t_max < Mesh.C(cell, 1)) && (Mesh.C(cell, 1) < 0.5 + t_max))
      phi_func(cell) = sin(3*pi*(Mesh.C(cell, 1)-1/6 - t_max))**2;
     endif
   endfor
endif


x = Mesh.C(:, 1);
plot (x, phi_k_minmod, "-o", x, phi_k_superbee, "-*", x, phi_k_upwind, "-+", x, phi_func, '-');

## Placing legend inside should return axes to original size
legend ({"Minmod", "SuperBee", "UpWind", "Exact"} , "location", "northwest");






