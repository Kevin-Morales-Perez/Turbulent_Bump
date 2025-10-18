function [vorticity] = computationVorticity(grad_u,grad_v,vorticity,ny,nx)
    %Computation of vorticity
    for i=1:ny
        for j=1:nx
            dv_dx=reshape(grad_v(i,j,1,1),[1,1]);
            du_dy=reshape(grad_u(i,j,1,2),[1,1]);
            vorticity(i,j)=dv_dx-du_dy;           
        end 
    end
end