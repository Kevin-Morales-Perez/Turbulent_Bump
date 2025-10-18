function [u_face,v_face] = face_vel_intRCnonrthgn2(u,v,p,u_face,v_face,...
    uVecNormFaces,distNeighbNods,uVecNeighbNods,weightDistFactors,...
    grad_p,aP,aPv,cellVols,nx_upstr,nx_dwnstr)
    %u: x  component velocity field , cell centered 
    %v: y component velocity field, cell centered
    %p: pressure cell centered
    %u_face: velocity normal to faces weast and east in the face centers
    %v_face:velocity normal to faces north and south in the face centers
    %VecCentNodFaceCents: Array with vectors from cell center to face
    %centers
    %uVecNormFaces:Array with vectors normal outward to each face of cells
    %distNeighbNods: Distances between cell centers and neigborhood nodes
    %uVecNeighbNods:Unitary vectors from cell center to nb cell centers
    %grad_u: Gradient for u 
    %grad_v: Gradient for v
    %aP:Coeffitient from momentum equations
    %aPv:Coeffitient from momentum equations for symetric Boundary
    %condition
    %cellVols: volumes of the cells
    %nx_upstr: number of cells in the free upwind zone
    %nx_dwnstr: number of cells in the free downstream zone 

    %Rie chow interpolation for face velocities 
    %   uf= 0.5*(uP + uA) + df*(Pp-Pa) - 0.5((VP/aP)*grad(p)p + 
    %  (VA/aA)*grad(p)A)
    %Orthogonal - and non orthogonal meshes
    
    %sizes
    size_field=size(u);
    nx=size_field(2);
    ny=size_field(1);
    
    % west and east face velocities  (u_face)
    
    for j=2:nx
        for i=1:ny
        
            %Compute average velocity_____________________________________
        
            %get vector from cell center to east face center
        
            %vec_p_fce =reshape(VecCentNodFaceCents(i,j-1,3,:),[1,2]);
        
            %get both velocity components at the cell center
        
            u_p=u(i,j-1);
        
            v_p=v(i,j-1);

            %group in a vector

            U_P=[u_p,v_p];

            %get both velocity components at the next cell center
        
            u_e=u(i,j);
        
            v_e=v(i,j);

            %group in a vector

            U_E=[u_e,v_e];
        
            %Get both gradients at the cell center 
        
            %grad_u_p = reshape(grad_u(i,j-1,:,:),[1,2]);
            %grad_v_p = reshape(grad_v(i,j-1,:,:),[1,2]);
        
            %get unitarsy vector normal to face east
        
            n_e=reshape(uVecNormFaces(i,j-1,3,:),[1,2]);
        
        
            %Average velocity equal to (up + grad(u)*v_p_fce,vp + grad(v)*v_p_fce)
            %*n

            %get weight factor at cell center for node E

            weight_E=reshape(weightDistFactors(i,j,3),[1,1]);
        
            u_face_av=weight_E*dot(U_P,n_e) + (1-weight_E)*dot(U_E,n_e);
        
            %use Rie Chow interpolation formula___________________________
            
            %get distance from cell center to east cell center
        
            dist_p_e =reshape(distNeighbNods(i,j-1,3,:),[1,1]);
        
            %get normal vector from cell center to east nb cell center 
        
            uvec_p_e=reshape(uVecNeighbNods(i,j-1,3,:),[1,2]);

            %Pressure gradients
            grad_p_p=reshape(grad_p(i,j-1,:,:),[1,2]);
            grad_p_e=reshape(grad_p(i,j,:,:),[1,2]);
        
            %Apply formula
            u_face(i,j)= u_face_av + 0.5*(cellVols(i,j-1)/aP(i,j-1) + ...
                cellVols(i,j)/aP(i,j))*((p(i,j-1)-p(i,j))/dist_p_e) - ...
                0.5*dot((cellVols(i,j-1)/aP(i,j-1))*grad_p_p + ...
                (cellVols(i,j)/aP(i,j))*grad_p_e,uvec_p_e);
       
        end
    end

    % North and South face velocities  (v_face)
    
    for i=2:ny-1
        for j=1:nx


            %Compute average velocity____________________________________
        
            %get vector from cell center to south face center
        
            %vec_p_fcs =reshape(VecCentNodFaceCents(i-1,j,4,:),[1,2]);
        
            %get both velocity components at the cell center
        
            u_p=u(i-1,j);
        
            v_p=v(i-1,j);

            %group in a vector

            U_P=[u_p,v_p];

            %get both velocity components at the cell center of the next
            %cell
        
            u_s=u(i,j);
        
            v_s=v(i,j);

            %group in a vector

            U_S=[u_s,v_s];
            
        
            %Get both gradients at the cell center 
        
            %grad_u_p = reshape(grad_u(i-1,j,:,:),[1,2]);
            %grad_v_p = reshape(grad_v(i-1,j,:,:),[1,2]);

            %get unitary vector normal to face south
            %negative sign to make positive the velocity
            n_s=-reshape(uVecNormFaces(i-1,j,4,:),[1,2]);

            %get weight factor for node s
            weight_S=reshape(weightDistFactors(i,j,4),[1,1]);

            %Average velocity equal to (up + grad(u)*v_p_fce,vp + grad(v)*v_p_fce)
            %*n

            v_face_av =weight_S*dot(U_P,n_s) + (1-weight_S)*dot(U_S,n_s);

        
            %get distance from cell center to south cell center
            dist_p_s = reshape(distNeighbNods(i-1,j,4,:), [1, 1]);
        
            %get normal vector from cell center to south cell center 
            uvec_p_s = reshape(uVecNeighbNods(i-1,j,4,:), [1, 2]);

            %Pressure gradients
            grad_p_p = reshape(grad_p(i-1,j,:,:), [1, 2]);
            grad_p_s = reshape(grad_p(i,j,:,:), [1, 2]);
        
            %Apply formula
            v_face(i,j) = v_face_av + 0.5*(cellVols(i-1,j)/aP(i-1,j) + ...
                cellVols(i,j)/aP(i,j))*((p(i-1,j)-p(i,j))/dist_p_s) - ...
                0.5*dot((cellVols(i-1,j)/aP(i-1,j))*grad_p_p + ...
                (cellVols(i,j)/aP(i,j))*grad_p_s, uvec_p_s);

        end
    end

    i=ny;

    for j=1:nx

        %Compute average velocity____________________________________
        
        %get vector from cell center to south face center
    
        %vec_p_fcs =reshape(VecCentNodFaceCents(i-1,j,4,:),[1,2]);
    
        %get both velocity components at the cell center
    
        u_p=u(i-1,j);
    
        v_p=v(i-1,j);

        %group in a vector

        U_P=[u_p,v_p];

        %get both velocity components at the cell center of the next
        %cell
    
        u_s=u(i,j);
    
        v_s=v(i,j);

        %group in a vector

        U_S=[u_s,v_s];
    
        %Get both gradients at the cell center 
    
        %grad_u_p = reshape(grad_u(i-1,j,:,:),[1,2]);
        %grad_v_p = reshape(grad_v(i-1,j,:,:),[1,2]);

        %get unitary vector normal to face south
        %negative sign to make positive the velocity
        n_s=-reshape(uVecNormFaces(i-1,j,4,:),[1,2]);

        %get weight factor for node s
        weight_S=reshape(weightDistFactors(i,j,4),[1,1]);

        %Average velocity equal to (up + grad(u)*v_p_fce,vp + grad(v)*v_p_fce)
        %*n

        v_face_av = weight_S*dot(U_P,n_s) + (1-weight_S)*dot(U_S,n_s);
    
        %get distance from cell center to south cell center
        dist_p_s = reshape(distNeighbNods(i-1,j,4,:), [1, 1]);
    
        %get normal vector from cell center to south cell center 
        uvec_p_s = reshape(uVecNeighbNods(i-1,j,4,:), [1, 2]);

        %Pressure gradients
        grad_p_p = reshape(grad_p(i-1,j,:,:), [1, 2]);
        grad_p_s = reshape(grad_p(i,j,:,:), [1, 2]);

        if (j<=nx_upstr) || j>(nx-nx_dwnstr)
            %use aPv coeffitients (Plate zone)
    
            %Apply formula
            v_face(i,j) = v_face_av + 0.5*(cellVols(i-1,j)/aP(i-1,j) + ...
                cellVols(i,j)/aPv(j))*((p(i-1,j)-p(i,j))/dist_p_s) - ...
                0.5*dot((cellVols(i-1,j)/aP(i-1,j))*grad_p_p + ...
                (cellVols(i,j)/aPv(j))*grad_p_s, uvec_p_s);
        else 
            %Apply formula
            v_face(i,j) = v_face_av + 0.5*(cellVols(i-1,j)/aP(i-1,j) + ...
                cellVols(i,j)/aP(i,j))*((p(i-1,j)-p(i,j))/dist_p_s) - ...
                0.5*dot((cellVols(i-1,j)/aP(i-1,j))*grad_p_p + ...
                (cellVols(i,j)/aP(i,j))*grad_p_s, uvec_p_s);

        end
        
    end
end