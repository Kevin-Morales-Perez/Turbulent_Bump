function [phiVert] = computePhiVertex(phi,phiGradients,...
    VecCentVertx,phiVert)
%Computation of transported variable phi
%  at vertex of each cell using gradients
%Can be used for Temperature , Velocity and less common Pressure 
%used for nu~ Spallart Allmaras transported variable 

    %Get the size of the mesh
    size_mesh=size(phi);
    nx=size_mesh(2);
    ny=size_mesh(1);

    for i=1:ny
        for j=1:nx

            %get gradient from cell i , j
            grad_phi=reshape(phiGradients(i,j,:,:),[1,2]);

            %get temperature from cell i, j

            phi_cell=phi(i,j);

            %get vectors from all vertex of cell i , j

            vecvert=reshape(VecCentVertx(i,j,:,:),[4,2]);

            %split vectors

            vec_wn=vecvert(1,:);
            vec_en=vecvert(2,:);
            vec_es=vecvert(3,:);
            vec_ws=vecvert(4,:);

            %compute temperatures at vertex

            phi_wn=phi_cell + dot(grad_phi,vec_wn);

            phi_en=phi_cell + dot(grad_phi,vec_en);

            phi_es=phi_cell + dot(grad_phi,vec_es);

            phi_ws=phi_cell + dot(grad_phi,vec_ws);

            %Store values 

            phiVert(i,j,:,:)=[phi_wn,phi_en,phi_es,phi_ws];

        end
    end

end
