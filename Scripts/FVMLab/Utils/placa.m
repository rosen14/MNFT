function [patches] = placa(xnod, faces, bf)
    % there are 4 patches, top, bottom, left, right
    % we look them geometrically

    [xmax] = max(xnod);
    [xmin] = min(xnod);

    f_left = abs(xnod(faces(bf,1),1)-xmin(1))<1e-3 & abs(xnod(faces(bf,2),1)-xmin(1))<1e-3; f_left = bf(f_left);
    f_right = abs(xnod(faces(bf,1),1)-xmax(1))<1e-3 & abs(xnod(faces(bf,2),1)-xmax(1))<1e-3; f_right = bf(f_right);
    f_top = abs(xnod(faces(bf,1),2)-xmax(2))<1e-3 & abs(xnod(faces(bf,2),2)-xmax(2))<1e-3; f_top = bf(f_top);
    f_bot = abs(xnod(faces(bf,1),2)-xmin(2))<1e-3 & abs(xnod(faces(bf,2),2)-xmin(2))<1e-3; f_bot = bf(f_bot);

<<<<<<< HEAD
    patches{1}.faces = f_left; patches{1}.name = 'left'; patches{1}.type = 'Dirichlet'; patches{1}.value = 1;
    patches{2}.faces = f_right; patches{2}.name = 'right'; patches{2}.type = 'Neumann'; patches{2}.value = 0;
    patches{3}.faces = f_top; patches{3}.name = 'top'; patches{3}.type = 'Neumann'; patches{3}.value = 0;
    patches{4}.faces = f_bot;patches{4}.name = 'bottom'; patches{4}.type = 'Dirichlet'; patches{4}.value = 0;
=======
    patches{1}.faces = f_left; patches{1}.name = 'left'; patches{1}.type = 'Dirichlet'; patches{1}.value = 2;
    patches{2}.faces = f_right; patches{2}.name = 'right'; patches{2}.type = 'Neumann'; patches{2}.value = 0;
    patches{3}.faces = f_top; patches{3}.name = 'top'; patches{3}.type = 'Neumann'; patches{3}.value = 0;
    patches{4}.faces = f_bot;patches{4}.name = 'bottom'; patches{4}.type = 'Dirichlet'; patches{4}.value = 1;
>>>>>>> 933214e (FVM HR)
end
