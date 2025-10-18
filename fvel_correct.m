function [u_face,v_face] = fvel_correct(u_face,v_face,p_prime,aP,aPv,...
    cellVols,weightDistFactors,distNeighbNods,nx_upstr,nx_dwnstr,alpha_uv)
%FACE VELOCITY CORRECTION
%u_face,v_face: velocities normal to faces
%p_prime: pressure correction field at cell centers
%%aP, aPv:central Coeffitient from momentum equations  
%cellVols: volumes of the cells
%weightDistFactors: Weight distance factors
%distNeighbNods: Distances between cell centers and neigborhood nodes
%nx_upstr: number of cells in the free upwind zone
%nx_dwnstr: number of cells in the free downstream zone

    %sizes
    size_field=size(p_prime);
    nx=size_field(2);
    ny=size_field(1);

    
    %Weast and East Faces

    for i = 1:ny
        for j=2:nx
            
            %Get weight distance factor to east node

            weight_e=reshape(weightDistFactors(i,j-1,3),[1,1]);

            %Get volume from cell P and cell E
            volP=cellVols(i,j-1);
            volE=cellVols(i,j);

            %get distance from node p to node e

            dist_e=reshape(distNeighbNods(i,j-1,3,:),[1,1]);

            u_face(i,j)=u_face(i,j) + alpha_uv*((volP/aP(i,j-1))*...
                weight_e + (1-weight_e)*(volE/aP(i,j)))*((...
                p_prime(i,j-1)-p_prime(i,j))/dist_e);
        end
    end

    %South and North faces 

    for i=2:ny
        for j=1:nx


            %Get weight distance factor to south node

            weight_s=reshape(weightDistFactors(i-1,j,4),[1,1]);

            %Get volume from cell P and cell S
            volP=cellVols(i-1,j);
            volS=cellVols(i,j);

            %get distance from node p to node s

            dist_s=reshape(distNeighbNods(i-1,j,4,:),[1,1]);
            
            if ((j<=nx_upstr) || j>(nx-nx_dwnstr)) && (i==ny)
                %Over the plate
                v_face(i,j)=v_face(i,j) +alpha_uv*((volP/aP(i-1,j))*...
                    weight_s + (1-weight_s)*(volS/aPv(j)))*...
                    ((p_prime(i-1,j)- p_prime(i,j))/dist_s);

            else

                v_face(i,j)=v_face(i,j) + alpha_uv*((volP/aP(i-1,j))*...
                    weight_s + (1-weight_s)*(volS/aP(i,j)))*...
                    ((p_prime(i-1,j)- p_prime(i,j))/dist_s);
            
            end
        end
    end
end