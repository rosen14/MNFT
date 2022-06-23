function phi_init = initial_condition(Mesh, choice)

   phi_init = zeros(Mesh.ncells,1);

   if strcmp(choice, 'suave')
     for cell = 1:Mesh.ncells
        if ((1/6 < Mesh.C(cell, 1)) && (Mesh.C(cell, 1) < 0.5))
          phi_init(cell) = sin(3*pi*(Mesh.C(cell, 1)-1/6))**2;
        endif
     endfor
   elseif strcmp(choice, 'step')
     for cell = 1:Mesh.ncells
        if ((1/6 < Mesh.C(cell, 1)) && (Mesh.C(cell, 1) < 0.5))
          phi_init(cell) = 1;
        endif
     endfor
   endif

end
