function [mu_turbulent_fw,mu_turbulent_fn,mu_turbulent_fe,...
    mu_turbulent_fs] = computePhiFaces(grad_mu_turbulent,mu_turbulent,...
    VecCentNodFaceCents,nx,ny,mu_turbulent_fw,mu_turbulent_fn,...
    mu_turbulent_fe,mu_turbulent_fs)
    %COMPUTE TRASNPORTED VARIABLE PHI AT FACES USING PHI AT CELL CENTROIDS 
    %AND GRADIENTS 
    
    
    for i=1:ny
        for j=1:nx
    
            %get gradient at cell center
    
            grad_phi=reshape(grad_mu_turbulent(i,j,:,:),[1,2]);
    
            %get temperature from cell i, j
    
            phi_cell=mu_turbulent(i,j);
    
    
            %get vectors from node to each center face;
    
            vecvert=reshape(VecCentNodFaceCents(i,j,:,:),[4,2]);
    
            %split vectors
    
            vec_wn=vecvert(1,:);
            vec_en=vecvert(2,:);
            vec_es=vecvert(3,:);
            vec_ws=vecvert(4,:);
    
    
    
            %face w_________________________________________________________
            mu_turbulent_fw(i,j)=phi_cell + dot(grad_phi,vec_wn);
    
            %face n_________________________________________________________
            mu_turbulent_fn(i,j)=phi_cell + dot(grad_phi,vec_en);
    
            %face e_________________________________________________________
            mu_turbulent_fe(i,j)=phi_cell + dot(grad_phi,vec_es);
    
            %face s_________________________________________________________
            mu_turbulent_fs(i,j)=phi_cell + dot(grad_phi,vec_ws);
    
        end 
    end


end