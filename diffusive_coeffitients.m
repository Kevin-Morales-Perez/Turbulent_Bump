function [dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c] = ...
    diffusive_coeffitients(nx,ny,distNeighbNods,uVecNormFaces,...
    uVecNeighbNods,lgtFaces,uVecParlFaces,mu,dW,dN,dE,dS,dW_c,dN_c,...
    dE_c,dS_c)
%COMPUTATION OF DIRECT AND NON ORTHOGONAL DIFFUSSION COEFFITIENTS 

%% DIRECT DIFFUSION
%face w difussive flux computation
for i=1:ny
    for j=1:nx
        %Get distance from node P to no node W
        dist_p_w=reshape(distNeighbNods(i,j,1),[1,1]);
        %get unitary vector normal to face w
        un_fw=reshape(uVecNormFaces(i,j,1,:),[1,2]);
        %get unitary vector from node p to node w
        uv_pw=reshape(uVecNeighbNods(i,j,1,:),[1,2]);
        %get length of face w
        lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);

        %Compute operator
        dW(i,j)=(mu*lgt_fw/(dist_p_w*dot(un_fw,uv_pw)));
        
    end
end


%face n difussive flux computation 
for i=1:ny
    for j=1:nx
        %Get distance from node P to no node N
        dist_p_n=reshape(distNeighbNods(i,j,2),[1,1]);
        %get unitary vector normal to face n
        un_fn=reshape(uVecNormFaces(i,j,2,:),[1,2]);
        %get unitary vector from node p to node n
        uv_pn=reshape(uVecNeighbNods(i,j,2,:),[1,2]);
        %get length of face n
        lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);

        %Compute operator
        dN(i,j)=(mu*lgt_fn/(dist_p_n*dot(un_fn,uv_pn)));
        
    end
end
%face e difussive flux computation
for i=1:ny
    for j=1:nx
        %Get distance from node P to no node e
        dist_p_e=reshape(distNeighbNods(i,j,3),[1,1]);
        %get unitary vector normal to face e
        un_fe=reshape(uVecNormFaces(i,j,3,:),[1,2]);
        %get unitary vector from node p to node e
        uv_pe=reshape(uVecNeighbNods(i,j,3,:),[1,2]);
        %get length of face e
        lgt_fe=reshape(lgtFaces(i,j,3),[1 1]);

        %Compute operator
        dE(i,j)=(mu*lgt_fe/(dist_p_e*dot(un_fe,uv_pe)));
        
    end
end

%face s difussive flux computation 
for i=1:ny
    for j=1:nx
        %Get distance from node P to no node S
        dist_p_s=reshape(distNeighbNods(i,j,4),[1,1]);
        %get unitary vector normal to face S
        un_fs=reshape(uVecNormFaces(i,j,4,:),[1,2]);
        %get unitary vector from node p to node s
        uv_ps=reshape(uVecNeighbNods(i,j,4,:),[1,2]);
        %get length of face s
        lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);

        %Compute operator
        dS(i,j)=(mu*lgt_fs/(dist_p_s*dot(un_fs,uv_ps)));
        
    end
end

%% NON ORTHOGONAL DIFFUSION

%face w
for i=1:ny
    for j=1:nx

        %get unitary vector parallel to face W
        up_fw=reshape(uVecParlFaces(i,j,1,:),[1,2]);
        %get unitary vector normal to face w
        un_fw=reshape(uVecNormFaces(i,j,1,:),[1,2]);
        %get unitary vector from node p to node w
        uv_pw=reshape(uVecNeighbNods(i,j,1,:),[1,2]);
        
        %Compute operator for cross diffusion terms
        dW_c(i,j) = -mu*(dot(up_fw,uv_pw)/dot(un_fw,uv_pw));
            
    end
end


%face n
for i=1:ny
    for j=1:nx

        %get unitary vector parallel to face N
        up_fn=reshape(uVecParlFaces(i,j,2,:),[1,2]);
        %get unitary vector normal to face n
        un_fn=reshape(uVecNormFaces(i,j,2,:),[1,2]);
        %get unitary vector from node p to node n
        uv_pn=reshape(uVecNeighbNods(i,j,2,:),[1,2]);
        

        %Compute operator for cross diffusion terms
        dN_c(i,j) = -mu*(dot(up_fn,uv_pn)/dot(un_fn,uv_pn));
            
    end
end

%face e
for i=1:ny
    for j=1:nx

        %get unitary vector parallel to face E
        up_fe=reshape(uVecParlFaces(i,j,3,:),[1,2]);
        %get unitary vector normal to face e
        un_fe=reshape(uVecNormFaces(i,j,3,:),[1,2]);
        %get unitary vector from node p to node e
        uv_pe=reshape(uVecNeighbNods(i,j,3,:),[1,2]);
        
        %Compute operator for cross diffusion terms
        dE_c(i,j) = -mu*(dot(up_fe,uv_pe)/dot(un_fe,uv_pe));
            
    end
end

%face s
for i=1:ny
    for j=1:nx

        %get unitary vector parallel to face S
        up_fs=reshape(uVecParlFaces(i,j,4,:),[1,2]);
        %get unitary vector normal to face s
        un_fs=reshape(uVecNormFaces(i,j,4,:),[1,2]);
        %get unitary vector from node p to node s
        uv_ps=reshape(uVecNeighbNods(i,j,4,:),[1,2]);
        
        %Compute operator for cross diffusion terms
        dS_c(i,j) = -mu*(dot(up_fs,uv_ps)/dot(un_fs,uv_ps));
            
    end
end

end