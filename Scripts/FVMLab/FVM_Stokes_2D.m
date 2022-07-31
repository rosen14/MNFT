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
BCs_u = @placa_u;
BCs_v = @placa_v;
BCs_p = @placa_p;
BCs_none = @placa_none;
% Difusividad
nu = 10; #Viscosidad cinemática m2/s
Q = 0;




%% END USER PARAMETERS %%%

% generacion de las estructuras geometricas
%make_data;
%filename = 'plano_triang_delaunay.msh';
filename = 'plano_regular.msh';
%filename = 'barra1D.msh';
Mesh = gen2DMesh(filename);


% patches para condiciones de borde
bf = find(Mesh.neighbour == 0);

patches_u = BCs_u(Mesh.xnod, Mesh.faces, bf);
patches_v = BCs_v(Mesh.xnod, Mesh.faces, bf);
patches_p = BCs_p(Mesh.xnod, Mesh.faces, bf);
patches_none = BCs_none(Mesh.xnod, Mesh.faces, bf); #Neumann con valor nulo

u_n = zeros(Mesh.ncells,1); #n = 0
v_n = zeros(Mesh.ncells,1);
P_k = zeros(Mesh.ncells,1);

[Ad_x, bd_x] = assemble_diffusion(Mesh, patches_u, nu);
[Ad_y, bd_y] = assemble_diffusion(Mesh, patches_v, nu);
A_x = Ad_x;
A_y = Ad_y;
wv = 0.7;
wp = 0.3;

d = eye(Mesh.ncells, 'logical');
nExtIt = 100;
nIntIt = 1;

A_x(d) = (1/wv).*A_x(d); %Verificar
A_y(d) = (1/wv).*A_y(d);

bs_x = assemble_source(Mesh, Q);
bs_y = assemble_source(Mesh, Q);
for iExtIt = 1:1:nExtIt

  iExtIt
  #Calculo gradiente de presión en celdas a partir de interpolar presión
  #en las caras y luego hacer gradiente Gauss
  [Pf, grad_P] = interp_phif_v2(Mesh, patches_p, P_k);



  bP_x = assemble_source(Mesh, - grad_P(:, 1));
  bP_y = assemble_source(Mesh, - grad_P(:, 2));

  b_relax_x = ((1-wv)/wv) * diag(Ad_x).*u_n;
  b_relax_y = ((1-wv)/wv) * diag(Ad_y).*v_n;



  b_x = b_relax_x + bs_x + bd_x + bP_x;
  b_y = b_relax_y + bs_y + bd_y + bP_y;


  u_n1 = A_x\b_x;
  v_n1 = A_y\b_y;
  for iIntit = 1:1:nIntIt
    iIntit
    #Matriz A promedio
    A = (A_x + A_y)/2;
    #Interpolar 1/aP en caras
    inv_aP = 1./(diag(A)./(Mesh.V));

    [inv_aPf, _] = interp_phif_v2(Mesh, patches_none, inv_aP);
    %
    %Computo H
    H_x = (-(A*u_n1 - diag(A).*u_n1) + bs_x + bd_x + b_relax_x)./(Mesh.V);
    H_y = (-(A*v_n1 - diag(A).*v_n1) + bs_y + bd_y + b_relax_y)./(Mesh.V);
    %Interpolar H/ap en caras (H/ap)f
    [inv_aP_Hf_x, _] = interp_phif_v2(Mesh, patches_none, inv_aP.*H_x);
    [inv_aP_Hf_y, _] = interp_phif_v2(Mesh, patches_none, inv_aP.*H_y);
    inv_aP_Hf = [inv_aP_Hf_x, inv_aP_Hf_y];
    [div_inv_aP_Hf] = 1.*fvc_div(inv_aP_Hf, Mesh); #Término derecho de la ec. de Presión
    %Ensamble de ecuación de presión. Análogo a difusión con fuente.
    [Ap, bp] = assemble_diffusion(Mesh, patches_p, -inv_aPf);
    bp = bp + div_inv_aP_Hf; #Ensamblo el lado derecho completo
    P_k1 = Ap\bp; #Presión en celdas

    %Corrección de flujos
    %[flux] = Eqn_flux(Mesh, patches, nu, A, phi)

    %Relajación de p
    P_k1 = wp*P_k1 + (1 - wp)*P_k;
    P_k = P_k1;
    [Pf, grad_P] = interp_phif_v2(Mesh, patches_p, P_k1);
    u_n1 = H_x.*inv_aP - grad_P(:, 1).*inv_aP;
    v_n1 = H_y.*inv_aP - grad_P(:, 2).*inv_aP;

    u_n = u_n1;
    v_n = v_n1;
  endfor

endfor

view2d_by_ele(Mesh.xnod, Mesh.icone, P_k);
axis equal;

colormap jet;
colorbar;

