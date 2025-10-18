function [u_star,v_star] = cvel_correct(aP,aPv,u,v,grad_p_prime,...
    cellVols,nx_upstr,nx_dwnstr,alpha_uv)
%CORRECTION OF CELL CENTER VELOCITIES 
% using pressure corection field to correct pressure velocities
%u_star, v_star: corrected velocities
%aP, aPv:Coeffitient from momentum equations  
%u,v: velocities at nodes
%grad_p_prime: gradient for pressure correction
%cellVols: volumes of the cells
%nx_upstr: number of cells in the free upwind zone
%nx_dwnstr: number of cells in the free downstream zone

    %sizes
    size_field=size(u);
    nx=size_field(2);
    ny=size_field(1);


    u_star=zeros(ny,nx);
    v_star=zeros(ny,nx);
    
    %Velocity in X and Y axis 

    for i = 1:ny-1
        for j = 1:nx
            
            %Get pressure correction gradient 
            gradP_p=reshape(grad_p_prime(i,j,:,:),[1,2]);

            gradP_px=gradP_p(1);
            gradP_py=gradP_p(2);
            
            u_star(i,j)=u(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_px;
 
            v_star(i,j)=v(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_py;
       
        end
    end

    i=ny;
    
    for j=1:nx
        %Get pressure correction gradient 
        gradP_p=reshape(grad_p_prime(i,j,:,:),[1,2]);

        gradP_px=gradP_p(1);
        gradP_py=gradP_p(2);

        u_star(i,j)=u(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_px;

        if (j<=nx_upstr) || j>(nx-nx_dwnstr)% Symetric BC  
            % coeffitients
            v_star(i,j)=v(i,j) - alpha_uv*(cellVols(i,j)/aPv(j))*gradP_py;
        else
            v_star(i,j)=v(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_py;
        end

    end

