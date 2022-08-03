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
nu = 0.03; #Viscosidad cinemática m2/s
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
flux = zeros(Mesh.nfaces, 1);

[Ad_x, bd_x] = assemble_diffusion(Mesh, patches_u, nu);
[Ad_y, bd_y] = assemble_diffusion(Mesh, patches_v, nu);

wv = 0.4;
wp = 0.6;

d = eye(Mesh.ncells, 'logical');
nExtIt = 100;
nIntIt = 1;


bs_x = assemble_source(Mesh, Q);
bs_y = assemble_source(Mesh, Q);
tvd_sch = 'upwind';#'central_diff';
for iExtIt = 1:1:nExtIt
  iExtIt

  [u_nf, grad_u_n] = interp_phif_v2(Mesh, patches_u, u_n);
  [v_nf, grad_v_n] = interp_phif_v2(Mesh, patches_v, v_n);
  vel = [u_nf, v_nf];
  [un_psif] = get_psif(Mesh, patches_u, u_n, vel, u_nf, grad_u_n, tvd_sch);
  [vn_psif] = get_psif(Mesh, patches_v, v_n, vel, v_nf, grad_v_n, tvd_sch);

  flux_w_vel = dot(vel(:, :),Mesh.Sf(:, :),2);
  #[Aa_x, ba_x] = assemble_advection_hr(Mesh, patches_u, vel, un_psif);
  #[Aa_y, ba_y] = assemble_advection_hr(Mesh, patches_v, vel, vn_psif);
  [Aa_x, ba_x] = assemble_advection_hr_flux(Mesh, patches_u, flux, un_psif);
  [Aa_y, ba_y] = assemble_advection_hr_flux(Mesh, patches_v, flux, vn_psif);

  A_x = Ad_x + 1.*Aa_x;
  A_y = Ad_y + 1.*Aa_y;

  b_relax_x = ((1-wv)/wv) * diag(A_x).*u_n; #Cuidado que este cálculo utiliza la matriz A sin relajar
  b_relax_y = ((1-wv)/wv) * diag(A_y).*v_n;

  A_x(d) = (1/wv).*A_x(d); #Cuidado que estoy relajando la variable de A, no hacer el calculo de b_relax luego de esto
  A_y(d) = (1/wv).*A_y(d);

  #Calculo gradiente de presión en celdas a partir de interpolar presión
  #en las caras y luego hacer gradiente Gauss
  [Pf, grad_P] = interp_phif_v2(Mesh, patches_p, P_k);

  bP_x = assemble_source(Mesh, - grad_P(:, 1));
  bP_y = assemble_source(Mesh, - grad_P(:, 2));

  b_x = b_relax_x + bs_x + bd_x + 1.*ba_x + bP_x ;
  b_y = b_relax_y + bs_y + bd_y + 1.*ba_y + bP_y;


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
    H_x = (-(A*u_n1 - diag(A).*u_n1) + bs_x + bd_x + 1.*ba_x + b_relax_x)./(Mesh.V);
    H_y = (-(A*v_n1 - diag(A).*v_n1) + bs_y + bd_y + 1.*ba_y + b_relax_y)./(Mesh.V);
    %Interpolar H/ap en caras (H/ap)f
    [inv_aP_Hf_x, _] = interp_phif_v2(Mesh, patches_none, inv_aP.*H_x);
    [inv_aP_Hf_y, _] = interp_phif_v2(Mesh, patches_none, inv_aP.*H_y);
    inv_aP_Hf = [inv_aP_Hf_x, inv_aP_Hf_y];
    [div_inv_aP_Hf] = 1.*fvc_div(inv_aP_Hf, Mesh); #Término derecho de la ec. de Presión
    %Ensamble de ecuación de presión. Análogo a difusión con fuente.
    [Ap, bp] = assemble_diffusion(Mesh, patches_p, -inv_aPf);
    #Fijo la presión en la celda 8 ya que el sistema es indeterminado (Presión
    #no fija en ninguna parte debido a Neumann en todas los bordes)

    bp = bp + div_inv_aP_Hf; #Ensamblo el lado derecho completo

    #bp(200) = bp(200) + Ap(200, 200)*0;
    #Ap(200, 200) = 2*Ap(200, 200);
    #El tema de fijar la presion en un punto creo que es la solucion al problema
    #Ap(150:150:end, 150:150:end) = Ap(150:150:end, 150:150:end) + 1000;

    P_k1 = Ap\bp; #Presión en celdas

    %Corrección de flujos
    [flux_p] = Eqn_flux(Mesh, patches_p, -inv_aPf, Ap, P_k1);
    flux = dot(inv_aP_Hf, Mesh.Sf, 2) - flux_p;
    [sum_flux_cell] = sum_flux_cells(flux, Mesh);
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

view2d_by_ele(Mesh.xnod, Mesh.icone, P_k1);
axis equal;

colormap jet;
colorbar;




#12 explota

