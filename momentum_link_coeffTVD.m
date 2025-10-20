function [aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeffTVD(nx,ny,...
    nx_upstr,nx_dwnstr,rho,lgtFaces,u_face,v_face,u0,...
    aW,aE,aN,aS,dW,dE,dN,dS,dW_c,dN_c,dE_c,dS_c,aP,aPv,suX,suY,...
    cellVols,u_corners,v_corners,grad_p,mu,mu_turbulent_fw,...
    mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs,grad_u,grad_v,...
    uVecNeighbNods,distNeighbNods,u,v)
%nx,ny: field sizes
%nx_upstr: number of cells in the free upwind zone
%nx_dwnstr: number of cells in the free downstream zone ss
%rho:density
%lgtFaces: Area of each face(lenghts in this case because is 2D)
%u_face: velocity normal to faces weast and east in the face centers
%v_face:velocity normal to faces north and south in the face centers
%u0: inlet velocity
%dW,dE,dN,dS:direct diffusion coeffitients
%dW_c,dN_c,dE_c,dS_c:Cross diffusion coeffitients
%aP,aPv:Central coeffitients
%suX,suY: Sources
%cellVols:Volume of each cell
%u_corners,v_corners: velocities at corners
%grad_p:Gradient of pressure 
%grad_tau_mn: Gradient of turbulent Stresses
%mu: dynamic viscosity
%mu_turbulent: Eddy viscosity


%Coeffitients for momentum equations 
%   Convection + Diffusion + Sources Coeffitiens for the momentum equations
%in order to have them in this format considering trasported
%  variable phi:
% Ap*phi_p = Aw*phi_w + An*phi_n + Ae*phi_e + As*phi_s + Su 

        %BOUNDARY  CONDITIONS

%       BOUNDARY CONDITIONS FOR FLOW OVER A BUMP
%           NEUMMAN  (DU/DY=0 ,DV/DY=0)
%                                                      O N  
%I                                                     U E   
%N                                                     T U  DU/DY=0
%L                                                     L M
%E                                                     E M  DV/DY=0
%T                                                     T A 
% ___ SYMETRIC ___|___   NO SLIP ___|___ SYMETRIC  ____  N

%nx_upstr,nx_fp,nx_dwnstr

    %_____________________________________________________________________
    %multipling each difussive term by dynamics viscosity plus eddy visc. 

    mu_total_fw=mu+ mu_turbulent_fw;
    mu_total_fn=mu+ mu_turbulent_fn;
    mu_total_fe=mu+ mu_turbulent_fe;
    mu_total_fs=mu+ mu_turbulent_fs;
    
    dW=dW.*(mu_total_fw);
    dE=dE.*(mu_total_fn);
    dN=dN.*(mu_total_fe);
    dS=dS.*(mu_total_fs);
    dW_c=dW_c.*(mu_total_fw);
    dN_c=dN_c.*(mu_total_fn);
    dE_c=dE_c.*(mu_total_fe);
    dS_c=dS_c.*(mu_total_fs);

    %_____________________________________________________________________


%1.- MOMENTUM LINK COEFFITIENTS
    
    %1.1- Interior Cells
    for i=2:ny-1
        for j=2:nx-1
            
            %get face areas
            lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
            lgt_fn=reshape(lgtFaces(i,j,2),[1 1]); 
            lgt_fe=reshape(lgtFaces(i,j,3),[1 1]); 
            lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);
    
            %Mass fluxes at faces 
            fW=-rho*lgt_fw*u_face(i,j);
            fE=rho*lgt_fe*u_face(i,j+1); 
            fN=rho*lgt_fn*v_face(i,j); 
            fS=-rho*lgt_fs*v_face(i+1,j); 
            
            %coeffitients using upwind squeme
            aW(i,j)=dW(i,j) + max(0,-fW);
            aE(i,j)=dE(i,j) + max(0,-fE);
            aN(i,j)=dN(i,j) + max(0,-fN);
            aS(i,j)=dS(i,j) + max(0,-fS);
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fW +fN +fE +fS;

            %explicit computation of cross - difussion components
            
            %______________________U__________________________
            
            %Get velocity vertexes values
            vert_u=reshape(u_corners(i,j,:,:),[1,4]);

            %Split velocities
            u_wn=vert_u(1);
            u_en=vert_u(2);
            u_es=vert_u(3);
            u_ws=vert_u(4);

            %compute cross diffusion terms for each face 
            Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
            Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
            Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
            Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
            
            %Sum all contributions
            Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

            %______________________V__________________________

            %Get velocity vertexes values
            vert_v=reshape(v_corners(i,j,:,:),[1,4]);

            %Split velocities
            v_wn=vert_v(1);
            v_en=vert_v(2);
            v_es=vert_v(3);
            v_ws=vert_v(4);

            %compute cross diffusion terms for each face 
            Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
            Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
            Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
            Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
            
            %Sum all contributions
            Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

                  
            
            %____________Get pressure Gradients__________________
            
            %get the gradient of pressure at the center of the cell
            gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
            
            
            %Grad in x
            grad_px=gradPcell(1);
            %Grad in y
            grad_py=gradPcell(2);

            %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

            %____________  TRANSPORTED VARIABLE = U   ______________________
                    
            %Get gradients at central node and neighborhod nodes 
            
            grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
            grad_u_W=reshape(grad_u(i,j-1,:,:),[1,2]);
            grad_u_N=reshape(grad_u(i-1,j,:,:),[1,2]);
            grad_u_E=reshape(grad_u(i,j+1,:,:),[1,2]);
            grad_u_S=reshape(grad_u(i+1,j,:,:),[1,2]);
            
             %____________  TRANSPORTED VARIABLE = V    _____________________
                        
            grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
            grad_v_W=reshape(grad_v(i,j-1,:,:),[1,2]);
            grad_v_N=reshape(grad_v(i-1,j,:,:),[1,2]);
            grad_v_E=reshape(grad_v(i,j+1,:,:),[1,2]);
            grad_v_S=reshape(grad_v(i+1,j,:,:),[1,2]);
            
            %vectors from central node_p to neighboorhod nodes (unitare vector
            %*magnitude)
            vec_P_W=reshape(uVecNeighbNods(i,j,1,:),[1,2])*distNeighbNods(i,j,1);
            vec_P_N=reshape(uVecNeighbNods(i,j,2,:),[1,2])*distNeighbNods(i,j,2);
            vec_P_E=reshape(uVecNeighbNods(i,j,3,:),[1,2])*distNeighbNods(i,j,3);
            vec_P_S=reshape(uVecNeighbNods(i,j,4,:),[1,2])*distNeighbNods(i,j,4);
            
            %face west********************************************************
            % mass flux fW
            
            %if possitive node  P is upwind and node W is downwind//////
            %if negative  node  W is upwind and node P is downwind////////
            
            if fW>0 %POSITIVE MASS FLUX
                %node  P upwind, node W downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_W;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i,j-1);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i,j-1);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_w=0;   
                else
                    r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_w_u=VanLeerGradientLimiter(r_w_u);
                
                    %compute deferred correction contribution for source
                    % in face W
                    suDC_w=-0.5*fW*psi_w_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_w=0;
                    
                else
                    r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_w_v=VanLeerGradientLimiter(r_w_v);
                
                    %compute deferred correction contribution for source
                    % in face W
                    svDC_w=-0.5*fW*psi_w_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  W upwind, node P downwind
                grad_phi_up_u=grad_u_W;
                grad_phi_up_v=grad_v_W;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_W;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j-1);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j-1);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_w=0;
                    
                else
                    r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_w_u=VanLeerGradientLimiter(r_w_u);
                    %compute deferred correction contribution for source
                    % in face W
                    suDC_w=-0.5*fW*psi_w_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_w=0;
                    
                else
                    r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_w_v=VanLeerGradientLimiter(r_w_v);
                    %compute deferred correction contribution for source
                    % in face W
                    svDC_w=-0.5*fW*psi_w_v*denom_v;
                end
            
            end
            
            %face north*******************************************************
            % mass flux fN
            %if possitive node  P is upwind and node N is downwind
            %if negative  node  N is upwind and node P is downwind
            
            if fN>0 %POSITIVE MASS FLUX
                %node  P upwind, node N downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_N;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i-1,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i-1,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_n=0;   
                else
                    r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_n_u=VanLeerGradientLimiter(r_n_u);
                
                    %compute deferred correction contribution for source
                    % in face N
                    suDC_n=-0.5*fN*psi_n_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_n=0;
                    
                else
                    r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_n_v=VanLeerGradientLimiter(r_n_v);
                
                    %compute deferred correction contribution for source
                    % in face N
                    svDC_n=-0.5*fN*psi_n_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  N upwind, node P downwind
                grad_phi_up_u=grad_u_N;
                grad_phi_up_v=grad_v_N;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_N;
                %upwind and downwind phi values for u
                phi_up_u=u(i-1,j);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i-1,j);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_n=0;
                    
                else
                    r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_n_u=VanLeerGradientLimiter(r_n_u);
                    %compute deferred correction contribution for source
                    % in face N
                    suDC_n=-0.5*fN*psi_n_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_n=0;
                    
                else
                    r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_n_v=VanLeerGradientLimiter(r_n_v);
                    %compute deferred correction contribution for source
                    % in face N
                    svDC_n=-0.5*fN*psi_n_v*denom_v;
                end
            
            end
            
            %face east********************************************************
            % mass flux fE
            %if possitive node P is upwind and node E is downwind
            %if negative  node E is upwind and node P is downwind
            
            if fE>0 %POSITIVE MASS FLUX
                %node  P upwind, node E downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_E;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i,j+1);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i,j+1);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_e=0;   
                else
                    r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_e_u=VanLeerGradientLimiter(r_e_u);
                
                    %compute deferred correction contribution for source
                    % in face E
                    suDC_e=-0.5*fE*psi_e_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_e=0;
                    
                else
                    r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_e_v=VanLeerGradientLimiter(r_e_v);
                
                    %compute deferred correction contribution for source
                    % in face E
                    svDC_e=-0.5*fE*psi_e_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  E upwind, node P downwind
                grad_phi_up_u=grad_u_E;
                grad_phi_up_v=grad_v_E;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_E;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j+1);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j+1);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_e=0;
                    
                else
                    r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_e_u=VanLeerGradientLimiter(r_e_u);
                    %compute deferred correction contribution for source
                    % in face E
                    suDC_e=-0.5*fE*psi_e_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_e=0;
                    
                else
                    r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_e_v=VanLeerGradientLimiter(r_e_v);
                    %compute deferred correction contribution for source
                    % in face E
                    svDC_e=-0.5*fE*psi_e_v*denom_v;
                end
            
            end
            
            %face south*******************************************************
            % mass flux fS
            %if possitive node  P is upwind and node S is downwind
            %if negative  node  S is upwind and node P is downwind
            
            if fS>0 %POSITIVE MASS FLUX
                %node  P upwind, node E downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_S;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i+1,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i+1,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_s=0;   
                else
                    r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_s_u=VanLeerGradientLimiter(r_s_u);
                
                    %compute deferred correction contribution for source
                    % in face S
                    suDC_s=-0.5*fS*psi_s_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_s=0;
                    
                else
                    r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_s_v=VanLeerGradientLimiter(r_s_v);
                
                    %compute deferred correction contribution for source
                    % in face S
                    svDC_s=-0.5*fS*psi_s_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  S upwind, node P downwind
                grad_phi_up_u=grad_u_S;
                grad_phi_up_v=grad_v_S;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_S;
                %upwind and downwind phi values for u
                phi_up_u=u(i+1,j);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i+1,j);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_s=0;
                    
                else
                    r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_s_u=VanLeerGradientLimiter(r_s_u);
                    %compute deferred correction contribution for source
                    % in face S
                    suDC_s=-0.5*fS*psi_s_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_s=0;
                    
                else
                    r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_s_v=VanLeerGradientLimiter(r_s_v);
                    %compute deferred correction contribution for source
                    % in face S
                    svDC_s=-0.5*fS*psi_s_v*denom_v;
                end
            
            end
            
            %total deferred correction
            %DC for U
            suDC_total= suDC_w + suDC_n + suDC_e + suDC_s;
            %DC for V
            svDC_total= svDC_w + svDC_n + svDC_e + svDC_s;

            %_____ Sources (Pressure derivative + cross difussion) ______

            suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + suDC_total;
            suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
    
        end
    end
    
    %1.2- Walls
    
    %1.2.1- Left wall - West (INLET)
    
    j=1;
    
    for i =2:ny-1

        %get face areas
        lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
        lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);
        lgt_fe=reshape(lgtFaces(i,j,3),[1 1]);
        lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);

        %Mass fluxes at faces 
        fW=-rho*lgt_fw*u_face(i,j); 
        fE=rho*lgt_fe*u_face(i,j+1); 
        fN=rho*lgt_fn*v_face(i,j); 
        fS=-rho*lgt_fs*v_face(i+1,j);   
            
        %coeffitients using upwind squeme
        aW(i,j)=dW(i,j); 
        aE(i,j)=dE(i,j) + max(0,-fE);
        aN(i,j)=dN(i,j) + max(0,-fN);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fN +fE +fS;

        %____________ CROSS DIFFUSION TERMS______________________
        %______________________U__________________________
            
        %Get velocity vertexes values
        vert_u=reshape(u_corners(i,j,:,:),[1,4]);

        %Split velocities
        u_wn=vert_u(1);
        u_en=vert_u(2);
        u_es=vert_u(3);
        u_ws=vert_u(4);

        %compute cross diffusion terms for each face 
        Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
        Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
        Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
        Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
        
        %Sum all contributions
        Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

            %______________________V__________________________

        %Get velocity vertexes values
        vert_v=reshape(v_corners(i,j,:,:),[1,4]);

        %Split velocities
        v_wn=vert_v(1);
        v_en=vert_v(2);
        v_es=vert_v(3);
        v_ws=vert_v(4);

        %compute cross diffusion terms for each face 
        Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
        Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
        Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
        Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
        
        %Sum all contributions
        Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

        %____________Get pressure Gradients__________________
            
        %get the gradient of pressure at the center of the cell
        gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
        
        %Grad in x
        grad_px=gradPcell(1);
        %Grad in y
        grad_py=gradPcell(2);

        %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

        %____________  TRANSPORTED VARIABLE = U   ______________________
                
        %Get gradients at central node and neighborhod nodes 
        
        grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
        grad_u_N=reshape(grad_u(i-1,j,:,:),[1,2]);
        grad_u_E=reshape(grad_u(i,j+1,:,:),[1,2]);
        grad_u_S=reshape(grad_u(i+1,j,:,:),[1,2]);
        
         %____________  TRANSPORTED VARIABLE = V    _____________________
                    
        grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
        grad_v_N=reshape(grad_v(i-1,j,:,:),[1,2]);
        grad_v_E=reshape(grad_v(i,j+1,:,:),[1,2]);
        grad_v_S=reshape(grad_v(i+1,j,:,:),[1,2]);
        
        %vectors from central node_p to neighboorhod nodes (unitare vector
        %*magnitude)
        vec_P_N=reshape(uVecNeighbNods(i,j,2,:),[1,2])*distNeighbNods(i,j,2);
        vec_P_E=reshape(uVecNeighbNods(i,j,3,:),[1,2])*distNeighbNods(i,j,3);
        vec_P_S=reshape(uVecNeighbNods(i,j,4,:),[1,2])*distNeighbNods(i,j,4);
        
        
        
        %face north*******************************************************
        % mass flux fN
        %if possitive node  P is upwind and node N is downwind
        %if negative  node  N is upwind and node P is downwind
        
        if fN>0 %POSITIVE MASS FLUX
            %node  P upwind, node N downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_N;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i-1,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i-1,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_n=0;   
            else
                r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_n_u=VanLeerGradientLimiter(r_n_u);
            
                %compute deferred correction contribution for source
                % in face N
                suDC_n=-0.5*fN*psi_n_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_n=0;
                
            else
                r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_n_v=VanLeerGradientLimiter(r_n_v);
            
                %compute deferred correction contribution for source
                % in face N
                svDC_n=-0.5*fN*psi_n_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  N upwind, node P downwind
            grad_phi_up_u=grad_u_N;
            grad_phi_up_v=grad_v_N;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_N;
            %upwind and downwind phi values for u
            phi_up_u=u(i-1,j);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i-1,j);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_n=0;
                
            else
                r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_n_u=VanLeerGradientLimiter(r_n_u);
                %compute deferred correction contribution for source
                % in face N
                suDC_n=-0.5*fN*psi_n_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_n=0;
                
            else
                r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_n_v=VanLeerGradientLimiter(r_n_v);
                %compute deferred correction contribution for source
                % in face N
                svDC_n=-0.5*fN*psi_n_v*denom_v;
            end
        
        end
        
        %face east********************************************************
        % mass flux fE
        %if possitive node P is upwind and node E is downwind
        %if negative  node E is upwind and node P is downwind
        
        if fE>0 %POSITIVE MASS FLUX
            %node  P upwind, node E downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_E;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i,j+1);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i,j+1);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_e=0;   
            else
                r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_e_u=VanLeerGradientLimiter(r_e_u);
            
                %compute deferred correction contribution for source
                % in face E
                suDC_e=-0.5*fE*psi_e_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_e=0;
                
            else
                r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_e_v=VanLeerGradientLimiter(r_e_v);
            
                %compute deferred correction contribution for source
                % in face E
                svDC_e=-0.5*fE*psi_e_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  E upwind, node P downwind
            grad_phi_up_u=grad_u_E;
            grad_phi_up_v=grad_v_E;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_E;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j+1);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j+1);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_e=0;
                
            else
                r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_e_u=VanLeerGradientLimiter(r_e_u);
                %compute deferred correction contribution for source
                % in face E
                suDC_e=-0.5*fE*psi_e_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_e=0;
                
            else
                r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_e_v=VanLeerGradientLimiter(r_e_v);
                %compute deferred correction contribution for source
                % in face E
                svDC_e=-0.5*fE*psi_e_v*denom_v;
            end
        
        end
        
        %face south*******************************************************
        % mass flux fS
        %if possitive node  P is upwind and node S is downwind
        %if negative  node  S is upwind and node P is downwind
        
        if fS>0 %POSITIVE MASS FLUX
            %node  P upwind, node E downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_S;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i+1,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i+1,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_s=0;   
            else
                r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_s_u=VanLeerGradientLimiter(r_s_u);
            
                %compute deferred correction contribution for source
                % in face S
                suDC_s=-0.5*fS*psi_s_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_s=0;
                
            else
                r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_s_v=VanLeerGradientLimiter(r_s_v);
            
                %compute deferred correction contribution for source
                % in face S
                svDC_s=-0.5*fS*psi_s_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  S upwind, node P downwind
            grad_phi_up_u=grad_u_S;
            grad_phi_up_v=grad_v_S;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_S;
            %upwind and downwind phi values for u
            phi_up_u=u(i+1,j);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i+1,j);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_s=0;
                
            else
                r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_s_u=VanLeerGradientLimiter(r_s_u);
                %compute deferred correction contribution for source
                % in face S
                suDC_s=-0.5*fS*psi_s_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_s=0;
                
            else
                r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_s_v=VanLeerGradientLimiter(r_s_v);
                %compute deferred correction contribution for source
                % in face S
                svDC_s=-0.5*fS*psi_s_v*denom_v;
            end
        
        end
        
        %total deferred correction
        %DC for U
        suDC_total=  + suDC_n + suDC_e + suDC_s;
        %DC for V
        svDC_total=  + svDC_n + svDC_e + svDC_s;

        %Sources (Pressure derivative)
        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(aW(i,j) -fW) + ...
            suDC_total;  
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
    
    end
    
    %1.2.2- Top wall - North  (NEUMMAN)
    
    i=1;
    
    for j=2:nx-1

         %get face areas
        lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
        lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);
        lgt_fe=reshape(lgtFaces(i,j,3),[1 1]);
        lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);
    
        %Mass fluxes
        fW=-rho*lgt_fw*u_face(i,j); 
        fE=rho*lgt_fe*u_face(i,j+1); 
        fN=rho*lgt_fn*v_face(i,j); 
        fS=-rho*lgt_fs*v_face(i+1,j);   
            
        %coeffitients using upwind scheme
        aW(i,j)=dW(i,j) + max(0,-fW);
        aE(i,j)=dE(i,j) + max(0,-fE);
        %aN(i,j)=dN(i,j);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aE(i,j)  + aS(i,j) + fW + fN + fE + fS;


        %explicit computation of cross - difussion components            
        %______________________U__________________________
        
        %Get velocity vertexes values
        vert_u=reshape(u_corners(i,j,:,:),[1,4]);

        %Split velocities
        u_wn=vert_u(1);
        u_en=vert_u(2);
        u_es=vert_u(3);
        u_ws=vert_u(4);

        %compute cross diffusion terms for each face 
        Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
        Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
        Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
        Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
        
        %Sum all contributions
        Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;
        
        %______________________v__________________________

        %Get velocity vertexes values
        vert_v=reshape(v_corners(i,j,:,:),[1,4]);

        %Split velocities
        v_wn=vert_v(1);
        v_en=vert_v(2);
        v_es=vert_v(3);
        v_ws=vert_v(4);

        %compute cross diffusion terms for each face 
        Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
        Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
        Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
        Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
        
        %Sum all contributions
        Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;
        
        
        
        %____________Compute pressure Gradients__________________

        %get the gradient of pressure at the center of the cell
        gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
        
        %Grad in x
        grad_px=gradPcell(1);
        %Grad in y
        grad_py=gradPcell(2);

        %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

        %____________  TRANSPORTED VARIABLE = U   ______________________
                
        %Get gradients at central node and neighborhod nodes 
        
        grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
        grad_u_W=reshape(grad_u(i,j-1,:,:),[1,2]);
        grad_u_E=reshape(grad_u(i,j+1,:,:),[1,2]);
        grad_u_S=reshape(grad_u(i+1,j,:,:),[1,2]);
        
         %____________  TRANSPORTED VARIABLE = V    _____________________
                    
        grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
        grad_v_W=reshape(grad_v(i,j-1,:,:),[1,2]);
        grad_v_E=reshape(grad_v(i,j+1,:,:),[1,2]);
        grad_v_S=reshape(grad_v(i+1,j,:,:),[1,2]);
        
        %vectors from central node_p to neighboorhod nodes (unitare vector
        %*magnitude)
        vec_P_W=reshape(uVecNeighbNods(i,j,1,:),[1,2])*distNeighbNods(i,j,1);
        vec_P_E=reshape(uVecNeighbNods(i,j,3,:),[1,2])*distNeighbNods(i,j,3);
        vec_P_S=reshape(uVecNeighbNods(i,j,4,:),[1,2])*distNeighbNods(i,j,4);
        
        %face west********************************************************
        % mass flux fW
        
        %if possitive node  P is upwind and node W is downwind//////
        %if negative  node  W is upwind and node P is downwind////////
        
        if fW>0 %POSITIVE MASS FLUX
            %node  P upwind, node W downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_W;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i,j-1);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i,j-1);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_w=0;   
            else
                r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_w_u=VanLeerGradientLimiter(r_w_u);
            
                %compute deferred correction contribution for source
                % in face W
                suDC_w=-0.5*fW*psi_w_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_w=0;
                
            else
                r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_w_v=VanLeerGradientLimiter(r_w_v);
            
                %compute deferred correction contribution for source
                % in face W
                svDC_w=-0.5*fW*psi_w_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  W upwind, node P downwind
            grad_phi_up_u=grad_u_W;
            grad_phi_up_v=grad_v_W;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_W;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j-1);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j-1);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_w=0;
                
            else
                r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_w_u=VanLeerGradientLimiter(r_w_u);
                %compute deferred correction contribution for source
                % in face W
                suDC_w=-0.5*fW*psi_w_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_w=0;
                
            else
                r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_w_v=VanLeerGradientLimiter(r_w_v);
                %compute deferred correction contribution for source
                % in face W
                svDC_w=-0.5*fW*psi_w_v*denom_v;
            end
        
        end
             
        %face east********************************************************
        % mass flux fE
        %if possitive node P is upwind and node E is downwind
        %if negative  node E is upwind and node P is downwind
        
        if fE>0 %POSITIVE MASS FLUX
            %node  P upwind, node E downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_E;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i,j+1);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i,j+1);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_e=0;   
            else
                r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_e_u=VanLeerGradientLimiter(r_e_u);
            
                %compute deferred correction contribution for source
                % in face E
                suDC_e=-0.5*fE*psi_e_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_e=0;
                
            else
                r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_e_v=VanLeerGradientLimiter(r_e_v);
            
                %compute deferred correction contribution for source
                % in face E
                svDC_e=-0.5*fE*psi_e_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  E upwind, node P downwind
            grad_phi_up_u=grad_u_E;
            grad_phi_up_v=grad_v_E;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_E;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j+1);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j+1);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_e=0;
                
            else
                r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_e_u=VanLeerGradientLimiter(r_e_u);
                %compute deferred correction contribution for source
                % in face E
                suDC_e=-0.5*fE*psi_e_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_e=0;
                
            else
                r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_e_v=VanLeerGradientLimiter(r_e_v);
                %compute deferred correction contribution for source
                % in face E
                svDC_e=-0.5*fE*psi_e_v*denom_v;
            end
        
        end
        
        %face south*******************************************************
        % mass flux fS
        %if possitive node  P is upwind and node S is downwind
        %if negative  node  S is upwind and node P is downwind
        
        if fS>0 %POSITIVE MASS FLUX
            %node  P upwind, node E downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_S;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i+1,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i+1,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_s=0;   
            else
                r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_s_u=VanLeerGradientLimiter(r_s_u);
            
                %compute deferred correction contribution for source
                % in face S
                suDC_s=-0.5*fS*psi_s_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_s=0;
                
            else
                r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_s_v=VanLeerGradientLimiter(r_s_v);
            
                %compute deferred correction contribution for source
                % in face S
                svDC_s=-0.5*fS*psi_s_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  S upwind, node P downwind
            grad_phi_up_u=grad_u_S;
            grad_phi_up_v=grad_v_S;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_S;
            %upwind and downwind phi values for u
            phi_up_u=u(i+1,j);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i+1,j);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_s=0;
                
            else
                r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_s_u=VanLeerGradientLimiter(r_s_u);
                %compute deferred correction contribution for source
                % in face S
                suDC_s=-0.5*fS*psi_s_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_s=0;
                
            else
                r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_s_v=VanLeerGradientLimiter(r_s_v);
                %compute deferred correction contribution for source
                % in face S
                svDC_s=-0.5*fS*psi_s_v*denom_v;
            end
        
        end
        
        %total deferred correction
        %DC for U
        suDC_total= suDC_w  + suDC_e + suDC_s;
        %DC for V
        svDC_total= svDC_w  + svDC_e + svDC_s;

        %Sources (Pressure dE(i,j)rivative + boundary condition)
        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*aN(i,j)+ suDC_total; 
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
    
    end
    
    %1.2.3- Right Wall - East  (Outlet)
    
    j=nx;
    for i=2:ny-1
       
        %get face areas
        lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
        lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);
        lgt_fe=reshape(lgtFaces(i,j,3),[1 1]);
        lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);

        %Mass fluxes at faces
        fW=-rho*lgt_fw*u_face(i,j);  
        fE=rho*lgt_fe*u_face(i,j+1); %Neumman condition 
        fN=rho*lgt_fn*v_face(i,j); 
        fS=-rho*lgt_fs*v_face(i+1,j);   
        
        %coeffitients using upwind squeme
        aW(i,j)=dW(i,j) + max(0,-fW);
        %aE(i,j)=0;
        aN(i,j)=dN(i,j) + max(0,-fN);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aN(i,j) + aS(i,j)  +fW +fN +fS +fE;

        %____________  CROSS DIFFUSION TERMS______________________

        %______________________U__________________________
            
        %Get velocity vertexes values
        vert_u=reshape(u_corners(i,j,:,:),[1,4]);

        %Split velocities
        u_wn=vert_u(1);
        u_en=vert_u(2);
        u_es=vert_u(3);
        u_ws=vert_u(4);

        %compute cross diffusion terms for each face 
        Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
        Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
        Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
        Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
        
        %Sum all contributions
        Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

        %______________________V__________________________

        %Get velocity vertexes values
        vert_v=reshape(v_corners(i,j,:,:),[1,4]);

        %Split velocities
        v_wn=vert_v(1);
        v_en=vert_v(2);
        v_es=vert_v(3);
        v_ws=vert_v(4);

        %compute cross diffusion terms for each face 
        Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
        Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
        Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
        Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
        
        %Sum all contributions
        Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

              
        
        %____________Get pressure Gradients__________________
        
        %get the gradient of pressure at the center of the cell
        gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
        
        
        %Grad in x
        grad_px=gradPcell(1);
        %Grad in y
        grad_py=gradPcell(2);

        %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

        %____________  TRANSPORTED VARIABLE = U   ______________________
                
        %Get gradients at central node and neighborhod nodes 
        
        grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
        grad_u_W=reshape(grad_u(i,j-1,:,:),[1,2]);
        grad_u_N=reshape(grad_u(i-1,j,:,:),[1,2]);
        grad_u_S=reshape(grad_u(i+1,j,:,:),[1,2]);
        
         %____________  TRANSPORTED VARIABLE = V    _____________________
                    
        grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
        grad_v_W=reshape(grad_v(i,j-1,:,:),[1,2]);
        grad_v_N=reshape(grad_v(i-1,j,:,:),[1,2]);
        grad_v_S=reshape(grad_v(i+1,j,:,:),[1,2]);
        
        %vectors from central node_p to neighboorhod nodes (unitare vector
        %*magnitude)
        vec_P_W=reshape(uVecNeighbNods(i,j,1,:),[1,2])*distNeighbNods(i,j,1);
        vec_P_N=reshape(uVecNeighbNods(i,j,2,:),[1,2])*distNeighbNods(i,j,2);
        vec_P_S=reshape(uVecNeighbNods(i,j,4,:),[1,2])*distNeighbNods(i,j,4);
        
        %face west********************************************************
        % mass flux fW
        
        %if possitive node  P is upwind and node W is downwind//////
        %if negative  node  W is upwind and node P is downwind////////
        
        if fW>0 %POSITIVE MASS FLUX
            %node  P upwind, node W downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_W;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i,j-1);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i,j-1);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_w=0;   
            else
                r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_w_u=VanLeerGradientLimiter(r_w_u);
            
                %compute deferred correction contribution for source
                % in face W
                suDC_w=-0.5*fW*psi_w_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_w=0;
                
            else
                r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_w_v=VanLeerGradientLimiter(r_w_v);
            
                %compute deferred correction contribution for source
                % in face W
                svDC_w=-0.5*fW*psi_w_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  W upwind, node P downwind
            grad_phi_up_u=grad_u_W;
            grad_phi_up_v=grad_v_W;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_W;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j-1);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j-1);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_w=0;
                
            else
                r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_w_u=VanLeerGradientLimiter(r_w_u);
                %compute deferred correction contribution for source
                % in face W
                suDC_w=-0.5*fW*psi_w_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_w=0;
                
            else
                r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_w_v=VanLeerGradientLimiter(r_w_v);
                %compute deferred correction contribution for source
                % in face W
                svDC_w=-0.5*fW*psi_w_v*denom_v;
            end
        
        end
        
        %face north*******************************************************
        % mass flux fN
        %if possitive node  P is upwind and node N is downwind
        %if negative  node  N is upwind and node P is downwind
        
        if fN>0 %POSITIVE MASS FLUX
            %node  P upwind, node N downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_N;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i-1,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i-1,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_n=0;   
            else
                r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_n_u=VanLeerGradientLimiter(r_n_u);
            
                %compute deferred correction contribution for source
                % in face N
                suDC_n=-0.5*fN*psi_n_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_n=0;
                
            else
                r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_n_v=VanLeerGradientLimiter(r_n_v);
            
                %compute deferred correction contribution for source
                % in face N
                svDC_n=-0.5*fN*psi_n_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  N upwind, node P downwind
            grad_phi_up_u=grad_u_N;
            grad_phi_up_v=grad_v_N;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_N;
            %upwind and downwind phi values for u
            phi_up_u=u(i-1,j);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i-1,j);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_n=0;
                
            else
                r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_n_u=VanLeerGradientLimiter(r_n_u);
                %compute deferred correction contribution for source
                % in face N
                suDC_n=-0.5*fN*psi_n_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_n=0;
                
            else
                r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_n_v=VanLeerGradientLimiter(r_n_v);
                %compute deferred correction contribution for source
                % in face N
                svDC_n=-0.5*fN*psi_n_v*denom_v;
            end
        
        end

        %face south*******************************************************
        % mass flux fS
        %if possitive node  P is upwind and node S is downwind
        %if negative  node  S is upwind and node P is downwind
        
        if fS>0 %POSITIVE MASS FLUX
            %node  P upwind, node E downwind
            grad_phi_up_u=grad_u_P;
            grad_phi_up_v=grad_v_P;
            %vector from upwind node to downwind node 
            up_do_vec=vec_P_S;
            %upwind and downwind phi values for u
            phi_up_u=u(i,j);
            phi_do_u=u(i+1,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i,j);
            phi_do_v=v(i+1,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
            
            %DEFERRED CORRECTION FOR U
            if denom_u==0 %DENOMINATOR IS ZERO
                suDC_s=0;   
            else
                r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_s_u=VanLeerGradientLimiter(r_s_u);
            
                %compute deferred correction contribution for source
                % in face S
                suDC_s=-0.5*fS*psi_s_u*denom_u;
                
            end
            %DEFERRED CORRECTION FOR V
            if denom_v==0
                svDC_s=0;
                
            else
                r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_s_v=VanLeerGradientLimiter(r_s_v);
            
                %compute deferred correction contribution for source
                % in face S
                svDC_s=-0.5*fS*psi_s_v*denom_v;
            end
        
        else %NEGATIVE MASS FLUX
            %node  S upwind, node P downwind
            grad_phi_up_u=grad_u_S;
            grad_phi_up_v=grad_v_S;
            %vector from upwind node to downwind node 
            up_do_vec=-vec_P_S;
            %upwind and downwind phi values for u
            phi_up_u=u(i+1,j);
            phi_do_u=u(i,j);
            %upwind and downwind phi values for v
            phi_up_v=v(i+1,j);
            phi_do_v=v(i,j);
            %denominators
            denom_u=phi_do_u -phi_up_u;
            denom_v=phi_do_v -phi_up_v;
        
            %DEFERRED CORRECTION FOR U
            if denom_u==0%DENOMINATOR IS ZERO
                suDC_s=0;
                
            else
                r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                psi_s_u=VanLeerGradientLimiter(r_s_u);
                %compute deferred correction contribution for source
                % in face S
                suDC_s=-0.5*fS*psi_s_u*denom_u;
            end
            
            %DEFERRED CORRECTION FOR v
            if denom_v==0%DENOMINATOR IS ZERO
                svDC_s=0;
                
            else
                r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                psi_s_v=VanLeerGradientLimiter(r_s_v);
                %compute deferred correction contribution for source
                % in face S
                svDC_s=-0.5*fS*psi_s_v*denom_v;
            end
        
        end
        
        %total deferred correction
        %DC for U
        suDC_total= suDC_w + suDC_n  + suDC_s;
        %DC for V
        svDC_total= svDC_w + svDC_n  + svDC_s;

        %_____ Sources (Pressure derivative + cross difussion) ______

        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + suDC_total;
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
    
    end
    
    %1.2.4- Bottom wall -- South (Symetric in the free-stream and no slip 
    % on the plate)
    
    i=ny;
    
    for j=2:nx-1

        %get face areas
        lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
        lgt_fn=reshape(lgtFaces(i,j,2),[1 1]); 
        lgt_fe=reshape(lgtFaces(i,j,3),[1 1]); 
        %lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);

        if (j<=nx_upstr) || j>(nx-nx_dwnstr)

            %central coeffitient for momentum equation for v will be
            %different due to Symetric BC (free- stream zones)

            %Mass fluxes at faces 
            fW=-rho*lgt_fw*u_face(i,j);
            fE=rho*lgt_fe*u_face(i,j+1); 
            fN=rho*lgt_fn*v_face(i,j); 
            %fS=-rho*lgt_fs*v_face(i+1,j); No mass flux at south face

            %coeffitients using upwind squeme
            aW(i,j)=dW(i,j) + max(0,-fW);
            aE(i,j)=dE(i,j) + max(0,-fE);
            aN(i,j)=dN(i,j) + max(0,-fN);
            aS(i,j)=dS(i,j);%No appllicable for U
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) +fW +fN +fE;
            aPv(j)=aW(i,j) + aE(i,j) + aN(i,j)+ aS(i,j) +fW +fN +fE;

            %____________ CROSS DIFFUSION TERMS______________________
            %______________________U__________________________
            
            %Get velocity vertexes values
            vert_u=reshape(u_corners(i,j,:,:),[1,4]);

            %Split velocities
            u_wn=vert_u(1);
            u_en=vert_u(2);
            u_es=vert_u(3);
            u_ws=vert_u(4);

            %compute cross diffusion terms for each face 
            Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
            Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
            Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
            Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
            
            %Sum all contributions
            Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

            %______________________V__________________________

            %Get velocity vertexes values
            vert_v=reshape(v_corners(i,j,:,:),[1,4]);

            %Split velocities
            v_wn=vert_v(1);
            v_en=vert_v(2);
            v_es=vert_v(3);
            v_ws=vert_v(4);

            %compute cross diffusion terms for each face 
            Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
            Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
            Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
            Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
            
            %Sum all contributions
            Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

            
            %____________Get pressure Gradients__________________
            
            %get the gradient of pressure at the center of the cell
            gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
            
            
            %Grad in x
            grad_px=gradPcell(1);
            %Grad in y
            grad_py=gradPcell(2);

            %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

            %____________  TRANSPORTED VARIABLE = U   ______________________
                    
            %Get gradients at central node and neighborhod nodes 
            
            grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
            grad_u_W=reshape(grad_u(i,j-1,:,:),[1,2]);
            grad_u_N=reshape(grad_u(i-1,j,:,:),[1,2]);
            grad_u_E=reshape(grad_u(i,j+1,:,:),[1,2]);
            
             %____________  TRANSPORTED VARIABLE = V    _____________________
                        
            grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
            grad_v_W=reshape(grad_v(i,j-1,:,:),[1,2]);
            grad_v_N=reshape(grad_v(i-1,j,:,:),[1,2]);
            grad_v_E=reshape(grad_v(i,j+1,:,:),[1,2]);
            
            %vectors from central node_p to neighboorhod nodes (unitare vector
            %*magnitude)
            vec_P_W=reshape(uVecNeighbNods(i,j,1,:),[1,2])*distNeighbNods(i,j,1);
            vec_P_N=reshape(uVecNeighbNods(i,j,2,:),[1,2])*distNeighbNods(i,j,2);
            vec_P_E=reshape(uVecNeighbNods(i,j,3,:),[1,2])*distNeighbNods(i,j,3);
            
            %face west********************************************************
            % mass flux fW
            
            %if possitive node  P is upwind and node W is downwind//////
            %if negative  node  W is upwind and node P is downwind////////
            
            if fW>0 %POSITIVE MASS FLUX
                %node  P upwind, node W downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_W;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i,j-1);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i,j-1);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_w=0;   
                else
                    r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_w_u=VanLeerGradientLimiter(r_w_u);
                
                    %compute deferred correction contribution for source
                    % in face W
                    suDC_w=-0.5*fW*psi_w_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_w=0;
                    
                else
                    r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_w_v=VanLeerGradientLimiter(r_w_v);
                
                    %compute deferred correction contribution for source
                    % in face W
                    svDC_w=-0.5*fW*psi_w_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  W upwind, node P downwind
                grad_phi_up_u=grad_u_W;
                grad_phi_up_v=grad_v_W;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_W;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j-1);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j-1);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_w=0;
                    
                else
                    r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_w_u=VanLeerGradientLimiter(r_w_u);
                    %compute deferred correction contribution for source
                    % in face W
                    suDC_w=-0.5*fW*psi_w_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_w=0;
                    
                else
                    r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_w_v=VanLeerGradientLimiter(r_w_v);
                    %compute deferred correction contribution for source
                    % in face W
                    svDC_w=-0.5*fW*psi_w_v*denom_v;
                end
            
            end
            
            %face north*******************************************************
            % mass flux fN
            %if possitive node  P is upwind and node N is downwind
            %if negative  node  N is upwind and node P is downwind
            
            if fN>0 %POSITIVE MASS FLUX
                %node  P upwind, node N downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_N;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i-1,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i-1,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_n=0;   
                else
                    r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_n_u=VanLeerGradientLimiter(r_n_u);
                
                    %compute deferred correction contribution for source
                    % in face N
                    suDC_n=-0.5*fN*psi_n_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_n=0;
                    
                else
                    r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_n_v=VanLeerGradientLimiter(r_n_v);
                
                    %compute deferred correction contribution for source
                    % in face N
                    svDC_n=-0.5*fN*psi_n_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  N upwind, node P downwind
                grad_phi_up_u=grad_u_N;
                grad_phi_up_v=grad_v_N;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_N;
                %upwind and downwind phi values for u
                phi_up_u=u(i-1,j);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i-1,j);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_n=0;
                    
                else
                    r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_n_u=VanLeerGradientLimiter(r_n_u);
                    %compute deferred correction contribution for source
                    % in face N
                    suDC_n=-0.5*fN*psi_n_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_n=0;
                    
                else
                    r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_n_v=VanLeerGradientLimiter(r_n_v);
                    %compute deferred correction contribution for source
                    % in face N
                    svDC_n=-0.5*fN*psi_n_v*denom_v;
                end
            
            end
            
            %face east********************************************************
            % mass flux fE
            %if possitive node P is upwind and node E is downwind
            %if negative  node E is upwind and node P is downwind
            
            if fE>0 %POSITIVE MASS FLUX
                %node  P upwind, node E downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_E;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i,j+1);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i,j+1);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_e=0;   
                else
                    r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_e_u=VanLeerGradientLimiter(r_e_u);
                
                    %compute deferred correction contribution for source
                    % in face E
                    suDC_e=-0.5*fE*psi_e_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_e=0;
                    
                else
                    r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_e_v=VanLeerGradientLimiter(r_e_v);
                
                    %compute deferred correction contribution for source
                    % in face E
                    svDC_e=-0.5*fE*psi_e_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  E upwind, node P downwind
                grad_phi_up_u=grad_u_E;
                grad_phi_up_v=grad_v_E;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_E;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j+1);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j+1);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_e=0;
                    
                else
                    r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_e_u=VanLeerGradientLimiter(r_e_u);
                    %compute deferred correction contribution for source
                    % in face E
                    suDC_e=-0.5*fE*psi_e_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_e=0;
                    
                else
                    r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_e_v=VanLeerGradientLimiter(r_e_v);
                    %compute deferred correction contribution for source
                    % in face E
                    svDC_e=-0.5*fE*psi_e_v*denom_v;
                end
            
            end
            
            %total deferred correction
            %DC for U
            suDC_total= suDC_w + suDC_n + suDC_e;
            %DC for V
            svDC_total= svDC_w + svDC_n + svDC_e;

          
            %_____ Sources (Pressure derivative + cross difussion) ______

            suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + suDC_total;
            suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;

        else

            %central Coeffitient for u and v are the same (over the plate)
    
            %Mass fluxes at faces 
            fW=-rho*lgt_fw*u_face(i,j);
            fE=rho*lgt_fe*u_face(i,j+1); 
            fN=rho*lgt_fn*v_face(i,j); 
            %fS=-rho*lgt_fs*v_face(i+1,j); No mass flux at south face
            
            %coeffitients using upwind squeme
            aW(i,j)=dW(i,j) + max(0,-fW);
            aE(i,j)=dE(i,j) + max(0,-fE);
            aN(i,j)=dN(i,j) + max(0,-fN);
            aS(i,j)=dS(i,j);
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fW +fN +fE;

            %explicit computation of cross - difussion components
            %______________________U__________________________
            
            %Get velocity vertexes values
            vert_u=reshape(u_corners(i,j,:,:),[1,4]);

            %Split velocities
            u_wn=vert_u(1);
            u_en=vert_u(2);
            u_es=vert_u(3);
            u_ws=vert_u(4);

            %compute cross diffusion terms for each face 
            Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
            Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
            Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
            Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
            
            %Sum all contributions
            Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

            %______________________V__________________________

            %Get velocity vertexes values
            vert_v=reshape(v_corners(i,j,:,:),[1,4]);

            %Split velocities
            v_wn=vert_v(1);
            v_en=vert_v(2);
            v_es=vert_v(3);
            v_ws=vert_v(4);

            %compute cross diffusion terms for each face 
            Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
            Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
            Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
            Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
            
            %Sum all contributions
            Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

                  
            
            %____________Get pressure Gradients__________________
            
            %get the gradient of pressure at the center of the cell
            gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
            
            
            %Grad in x
            grad_px=gradPcell(1);
            %Grad in y
            grad_py=gradPcell(2);

            %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

            %____________  TRANSPORTED VARIABLE = U   ______________________
                    
            %Get gradients at central node and neighborhod nodes 
            
            grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
            grad_u_W=reshape(grad_u(i,j-1,:,:),[1,2]);
            grad_u_N=reshape(grad_u(i-1,j,:,:),[1,2]);
            grad_u_E=reshape(grad_u(i,j+1,:,:),[1,2]);
            
             %____________  TRANSPORTED VARIABLE = V    _____________________
                        
            grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
            grad_v_W=reshape(grad_v(i,j-1,:,:),[1,2]);
            grad_v_N=reshape(grad_v(i-1,j,:,:),[1,2]);
            grad_v_E=reshape(grad_v(i,j+1,:,:),[1,2]);
            
            %vectors from central node_p to neighboorhod nodes (unitare vector
            %*magnitude)
            vec_P_W=reshape(uVecNeighbNods(i,j,1,:),[1,2])*distNeighbNods(i,j,1);
            vec_P_N=reshape(uVecNeighbNods(i,j,2,:),[1,2])*distNeighbNods(i,j,2);
            vec_P_E=reshape(uVecNeighbNods(i,j,3,:),[1,2])*distNeighbNods(i,j,3);
            
            %face west********************************************************
            % mass flux fW
            
            %if possitive node  P is upwind and node W is downwind//////
            %if negative  node  W is upwind and node P is downwind////////
            
            if fW>0 %POSITIVE MASS FLUX
                %node  P upwind, node W downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_W;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i,j-1);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i,j-1);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_w=0;   
                else
                    r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_w_u=VanLeerGradientLimiter(r_w_u);
                
                    %compute deferred correction contribution for source
                    % in face W
                    suDC_w=-0.5*fW*psi_w_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_w=0;
                    
                else
                    r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_w_v=VanLeerGradientLimiter(r_w_v);
                
                    %compute deferred correction contribution for source
                    % in face W
                    svDC_w=-0.5*fW*psi_w_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  W upwind, node P downwind
                grad_phi_up_u=grad_u_W;
                grad_phi_up_v=grad_v_W;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_W;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j-1);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j-1);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_w=0;
                    
                else
                    r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_w_u=VanLeerGradientLimiter(r_w_u);
                    %compute deferred correction contribution for source
                    % in face W
                    suDC_w=-0.5*fW*psi_w_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_w=0;
                    
                else
                    r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_w_v=VanLeerGradientLimiter(r_w_v);
                    %compute deferred correction contribution for source
                    % in face W
                    svDC_w=-0.5*fW*psi_w_v*denom_v;
                end
            
            end
            
            %face north*******************************************************
            % mass flux fN
            %if possitive node  P is upwind and node N is downwind
            %if negative  node  N is upwind and node P is downwind
            
            if fN>0 %POSITIVE MASS FLUX
                %node  P upwind, node N downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_N;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i-1,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i-1,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_n=0;   
                else
                    r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_n_u=VanLeerGradientLimiter(r_n_u);
                
                    %compute deferred correction contribution for source
                    % in face N
                    suDC_n=-0.5*fN*psi_n_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_n=0;
                    
                else
                    r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_n_v=VanLeerGradientLimiter(r_n_v);
                
                    %compute deferred correction contribution for source
                    % in face N
                    svDC_n=-0.5*fN*psi_n_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  N upwind, node P downwind
                grad_phi_up_u=grad_u_N;
                grad_phi_up_v=grad_v_N;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_N;
                %upwind and downwind phi values for u
                phi_up_u=u(i-1,j);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i-1,j);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_n=0;
                    
                else
                    r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_n_u=VanLeerGradientLimiter(r_n_u);
                    %compute deferred correction contribution for source
                    % in face N
                    suDC_n=-0.5*fN*psi_n_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_n=0;
                    
                else
                    r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_n_v=VanLeerGradientLimiter(r_n_v);
                    %compute deferred correction contribution for source
                    % in face N
                    svDC_n=-0.5*fN*psi_n_v*denom_v;
                end
            
            end
            
            %face east********************************************************
            % mass flux fE
            %if possitive node P is upwind and node E is downwind
            %if negative  node E is upwind and node P is downwind
            
            if fE>0 %POSITIVE MASS FLUX
                %node  P upwind, node E downwind
                grad_phi_up_u=grad_u_P;
                grad_phi_up_v=grad_v_P;
                %vector from upwind node to downwind node 
                up_do_vec=vec_P_E;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j);
                phi_do_u=u(i,j+1);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j);
                phi_do_v=v(i,j+1);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
                
                %DEFERRED CORRECTION FOR U
                if denom_u==0 %DENOMINATOR IS ZERO
                    suDC_e=0;   
                else
                    r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_e_u=VanLeerGradientLimiter(r_e_u);
                
                    %compute deferred correction contribution for source
                    % in face E
                    suDC_e=-0.5*fE*psi_e_u*denom_u;
                    
                end
                %DEFERRED CORRECTION FOR V
                if denom_v==0
                    svDC_e=0;
                    
                else
                    r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_e_v=VanLeerGradientLimiter(r_e_v);
                
                    %compute deferred correction contribution for source
                    % in face E
                    svDC_e=-0.5*fE*psi_e_v*denom_v;
                end
            
            else %NEGATIVE MASS FLUX
                %node  E upwind, node P downwind
                grad_phi_up_u=grad_u_E;
                grad_phi_up_v=grad_v_E;
                %vector from upwind node to downwind node 
                up_do_vec=-vec_P_E;
                %upwind and downwind phi values for u
                phi_up_u=u(i,j+1);
                phi_do_u=u(i,j);
                %upwind and downwind phi values for v
                phi_up_v=v(i,j+1);
                phi_do_v=v(i,j);
                %denominators
                denom_u=phi_do_u -phi_up_u;
                denom_v=phi_do_v -phi_up_v;
            
                %DEFERRED CORRECTION FOR U
                if denom_u==0%DENOMINATOR IS ZERO
                    suDC_e=0;
                    
                else
                    r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
                    psi_e_u=VanLeerGradientLimiter(r_e_u);
                    %compute deferred correction contribution for source
                    % in face E
                    suDC_e=-0.5*fE*psi_e_u*denom_u;
                end
                
                %DEFERRED CORRECTION FOR v
                if denom_v==0%DENOMINATOR IS ZERO
                    svDC_e=0;
                    
                else
                    r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
                    psi_e_v=VanLeerGradientLimiter(r_e_v);
                    %compute deferred correction contribution for source
                    % in face E
                    svDC_e=-0.5*fE*psi_e_v*denom_v;
                end
            
            end
            
            %total deferred correction
            %DC for U
            suDC_total= suDC_w + suDC_n + suDC_e;
            %DC for V
            svDC_total= svDC_w + svDC_n + svDC_e;

            suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + suDC_total;
            suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;

        end
    
    end
    
    %1.3- Corners
    
    %1.3.1- North-West *************************
    i=1;
    j=1;

    %get face areas
    lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
    lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);  
    lgt_fe=reshape(lgtFaces(i,j,3),[1 1]); 
    lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);
    
    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_face(i,j); 
    fE=rho*lgt_fe*u_face(i,j+1); 
    fN=rho*lgt_fn*v_face(i,j); 
    fS=-rho*lgt_fs*v_face(i+1,j);
    
    %coeffitients using upwind squeme
    aW(i,j)=dW(i,j);
    aE(i,j)=dE(i,j) + max(0,-fE);
    %aN(i,j)=dN(i,j);
    aS(i,j)=dS(i,j) + max(0,-fS);
    aP(i,j)=aW(i,j) + aE(i,j) + aS(i,j) + fN + fE +fS;
    
    %____________ CROSS DIFFUSION TERMS______________________
    
    %______________________U__________________________
            
    %Get velocity vertexes values
    vert_u=reshape(u_corners(i,j,:,:),[1,4]);

    %Split velocities
    u_wn=vert_u(1);
    u_en=vert_u(2);
    u_es=vert_u(3);
    u_ws=vert_u(4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values
    vert_v=reshape(v_corners(i,j,:,:),[1,4]);

    %Split velocities
    v_wn=vert_v(1);
    v_en=vert_v(2);
    v_es=vert_v(3);
    v_ws=vert_v(4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell
    gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
    
    
    %Grad in x
    grad_px=gradPcell(1);
    %Grad in y
    grad_py=gradPcell(2);

    %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

    %____________  TRANSPORTED VARIABLE = U   ______________________
            
    %Get gradients at central node and neighborhod nodes 
    
    grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
    grad_u_E=reshape(grad_u(i,j+1,:,:),[1,2]);
    grad_u_S=reshape(grad_u(i+1,j,:,:),[1,2]);
    
     %____________  TRANSPORTED VARIABLE = V    _____________________
                
    grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
    grad_v_E=reshape(grad_v(i,j+1,:,:),[1,2]);
    grad_v_S=reshape(grad_v(i+1,j,:,:),[1,2]);
    
    %vectors from central node_p to neighboorhod nodes (unitare vector
    %*magnitude)
    vec_P_E=reshape(uVecNeighbNods(i,j,3,:),[1,2])*distNeighbNods(i,j,3);
    vec_P_S=reshape(uVecNeighbNods(i,j,4,:),[1,2])*distNeighbNods(i,j,4);
    
    
    
    %face east********************************************************
    % mass flux fE
    %if possitive node P is upwind and node E is downwind
    %if negative  node E is upwind and node P is downwind
    
    if fE>0 %POSITIVE MASS FLUX
        %node  P upwind, node E downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_E;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i,j+1);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i,j+1);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_e=0;   
        else
            r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_e_u=VanLeerGradientLimiter(r_e_u);
        
            %compute deferred correction contribution for source
            % in face E
            suDC_e=-0.5*fE*psi_e_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_e=0;
            
        else
            r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_e_v=VanLeerGradientLimiter(r_e_v);
        
            %compute deferred correction contribution for source
            % in face E
            svDC_e=-0.5*fE*psi_e_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  E upwind, node P downwind
        grad_phi_up_u=grad_u_E;
        grad_phi_up_v=grad_v_E;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_E;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j+1);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j+1);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_e=0;
            
        else
            r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_e_u=VanLeerGradientLimiter(r_e_u);
            %compute deferred correction contribution for source
            % in face E
            suDC_e=-0.5*fE*psi_e_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_e=0;
            
        else
            r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_e_v=VanLeerGradientLimiter(r_e_v);
            %compute deferred correction contribution for source
            % in face E
            svDC_e=-0.5*fE*psi_e_v*denom_v;
        end
    
    end
    
    %face south*******************************************************
    % mass flux fS
    %if possitive node  P is upwind and node S is downwind
    %if negative  node  S is upwind and node P is downwind
    
    if fS>0 %POSITIVE MASS FLUX
        %node  P upwind, node E downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_S;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i+1,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i+1,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_s=0;   
        else
            r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_s_u=VanLeerGradientLimiter(r_s_u);
        
            %compute deferred correction contribution for source
            % in face S
            suDC_s=-0.5*fS*psi_s_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_s=0;
            
        else
            r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_s_v=VanLeerGradientLimiter(r_s_v);
        
            %compute deferred correction contribution for source
            % in face S
            svDC_s=-0.5*fS*psi_s_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  S upwind, node P downwind
        grad_phi_up_u=grad_u_S;
        grad_phi_up_v=grad_v_S;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_S;
        %upwind and downwind phi values for u
        phi_up_u=u(i+1,j);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i+1,j);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_s=0;
            
        else
            r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_s_u=VanLeerGradientLimiter(r_s_u);
            %compute deferred correction contribution for source
            % in face S
            suDC_s=-0.5*fS*psi_s_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_s=0;
            
        else
            r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_s_v=VanLeerGradientLimiter(r_s_v);
            %compute deferred correction contribution for source
            % in face S
            svDC_s=-0.5*fS*psi_s_v*denom_v;
        end
    
    end
    
    %total deferred correction
    %DC for U
    suDC_total= suDC_e + suDC_s;
    %DC for V
    svDC_total= svDC_e + svDC_s;

    %Sources (Pressure dErivative + boundary condition)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(aW(i,j) + aN(i,j) -fW)...
        + suDC_total;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
    
    %1.3.2- North-East ************************
    i=1;
    j=nx;
    
    %get face areas
    lgt_fw=reshape(lgtFaces(i,j,1),[1 1]); 
    lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);
    lgt_fe=reshape(lgtFaces(i,j,3),[1 1]);  
    lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);

    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_face(i,j);
    fE=rho*lgt_fe*u_face(i,j+1);
    fN=rho*lgt_fn*v_face(i,j); 
    fS=-rho*lgt_fs*v_face(i+1,j); 
    
    %coeffitients using upwind squeme
    aW(i,j)=dW(i,j) + max(0,-fW);
    %aE(i,j)=2*dE(i,j);
    %aN(i,j)=dN(i,j);
    aS(i,j)=dS(i,j) + max(0,-fS);
    aP(i,j)=aW(i,j) + aN(i,j) + aS(i,j) +fW + fN+fS +fE;

    %____________ CROSS DIFFUSION TERMS______________________

    %______________________U__________________________
            
    %Get velocity vertexes values
    vert_u=reshape(u_corners(i,j,:,:),[1,4]);

    %Split velocities
    u_wn=vert_u(1);
    u_en=vert_u(2);
    u_es=vert_u(3);
    u_ws=vert_u(4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values
    vert_v=reshape(v_corners(i,j,:,:),[1,4]);

    %Split velocities
    v_wn=vert_v(1);
    v_en=vert_v(2);
    v_es=vert_v(3);
    v_ws=vert_v(4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell
    gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
    
    
    %Grad in x
    grad_px=gradPcell(1);
    %Grad in y
    grad_py=gradPcell(2);

    %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

    %____________  TRANSPORTED VARIABLE = U   ______________________
            
    %Get gradients at central node and neighborhod nodes 
    
    grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
    grad_u_W=reshape(grad_u(i,j-1,:,:),[1,2]);
    grad_u_S=reshape(grad_u(i+1,j,:,:),[1,2]);
    
     %____________  TRANSPORTED VARIABLE = V    _____________________
                
    grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
    grad_v_W=reshape(grad_v(i,j-1,:,:),[1,2]);
    grad_v_S=reshape(grad_v(i+1,j,:,:),[1,2]);
    
    %vectors from central node_p to neighboorhod nodes (unitare vector
    %*magnitude)
    vec_P_W=reshape(uVecNeighbNods(i,j,1,:),[1,2])*distNeighbNods(i,j,1);
    vec_P_S=reshape(uVecNeighbNods(i,j,4,:),[1,2])*distNeighbNods(i,j,4);
    
    %face west********************************************************
    % mass flux fW
    
    %if possitive node  P is upwind and node W is downwind//////
    %if negative  node  W is upwind and node P is downwind////////
    
    if fW>0 %POSITIVE MASS FLUX
        %node  P upwind, node W downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_W;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i,j-1);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i,j-1);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_w=0;   
        else
            r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_w_u=VanLeerGradientLimiter(r_w_u);
        
            %compute deferred correction contribution for source
            % in face W
            suDC_w=-0.5*fW*psi_w_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_w=0;
            
        else
            r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_w_v=VanLeerGradientLimiter(r_w_v);
        
            %compute deferred correction contribution for source
            % in face W
            svDC_w=-0.5*fW*psi_w_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  W upwind, node P downwind
        grad_phi_up_u=grad_u_W;
        grad_phi_up_v=grad_v_W;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_W;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j-1);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j-1);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_w=0;
            
        else
            r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_w_u=VanLeerGradientLimiter(r_w_u);
            %compute deferred correction contribution for source
            % in face W
            suDC_w=-0.5*fW*psi_w_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_w=0;
            
        else
            r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_w_v=VanLeerGradientLimiter(r_w_v);
            %compute deferred correction contribution for source
            % in face W
            svDC_w=-0.5*fW*psi_w_v*denom_v;
        end
    
    end
    
    %face south*******************************************************
    % mass flux fS
    %if possitive node  P is upwind and node S is downwind
    %if negative  node  S is upwind and node P is downwind
    
    if fS>0 %POSITIVE MASS FLUX
        %node  P upwind, node E downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_S;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i+1,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i+1,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_s=0;   
        else
            r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_s_u=VanLeerGradientLimiter(r_s_u);
        
            %compute deferred correction contribution for source
            % in face S
            suDC_s=-0.5*fS*psi_s_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_s=0;
            
        else
            r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_s_v=VanLeerGradientLimiter(r_s_v);
        
            %compute deferred correction contribution for source
            % in face S
            svDC_s=-0.5*fS*psi_s_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  S upwind, node P downwind
        grad_phi_up_u=grad_u_S;
        grad_phi_up_v=grad_v_S;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_S;
        %upwind and downwind phi values for u
        phi_up_u=u(i+1,j);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i+1,j);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_s=0;
            
        else
            r_s_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_s_u=VanLeerGradientLimiter(r_s_u);
            %compute deferred correction contribution for source
            % in face S
            suDC_s=-0.5*fS*psi_s_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_s=0;
            
        else
            r_s_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_s_v=VanLeerGradientLimiter(r_s_v);
            %compute deferred correction contribution for source
            % in face S
            svDC_s=-0.5*fS*psi_s_v*denom_v;
        end
    
    end
    
    %total deferred correction
    %DC for U
    suDC_total= suDC_w  + suDC_s;
    %DC for V
    svDC_total= svDC_w  + svDC_s;

    %Sources (Pressure dE(i,j)rivative + boundary condition)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*aN(i,j) + suDC_total;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
    
    %1.3.3- South-East (Symetric at South and outlet at East)***********
    i=ny;
    j=nx;

    %get face areas
    lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
    lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);
    lgt_fe=reshape(lgtFaces(i,j,3),[1 1]);  
    %lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);
    
    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_face(i,j);  
    fE=rho*lgt_fe*u_face(i,j+1); 
    fN=rho*lgt_fn*v_face(i,j); 
    %fS=-rho*lgt_fs*v_face(i+1,j); No mass flux at south face
    
    %coeffitients using upwind scheme
    aW(i,j)=dW(i,j) + max(0,-fW);
    %aE(i,j)=2*dE(i,j);
    aN(i,j)=dN(i,j) + max(0,-fN);
    aS(i,j)=dS(i,j);%No applicable for u due to SBC
    aP(i,j)=aW(i,j) + aN(i,j) + fW + fN + fE;
    aPv(j)= aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fW +fN +fE;

    %____________ CROSS DIFFUSION TERMS______________________

    %______________________U__________________________
            
    %Get velocity vertexes values
    vert_u=reshape(u_corners(i,j,:,:),[1,4]);

    %Split velocities
    u_wn=vert_u(1);
    u_en=vert_u(2);
    u_es=vert_u(3);
    u_ws=vert_u(4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values
    vert_v=reshape(v_corners(i,j,:,:),[1,4]);

    %Split velocities
    v_wn=vert_v(1);
    v_en=vert_v(2);
    v_es=vert_v(3);
    v_ws=vert_v(4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell
    gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
    
    %Grad in x
    grad_px=gradPcell(1);
    %Grad in y
    grad_py=gradPcell(2);

    %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

    %____________  TRANSPORTED VARIABLE = U   ______________________
            
    %Get gradients at central node and neighborhod nodes 
    
    grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
    grad_u_W=reshape(grad_u(i,j-1,:,:),[1,2]);
    grad_u_N=reshape(grad_u(i-1,j,:,:),[1,2]);
    
     %____________  TRANSPORTED VARIABLE = V    _____________________
                
    grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
    grad_v_W=reshape(grad_v(i,j-1,:,:),[1,2]);
    grad_v_N=reshape(grad_v(i-1,j,:,:),[1,2]);
    
    %vectors from central node_p to neighboorhod nodes (unitare vector
    %*magnitude)
    vec_P_W=reshape(uVecNeighbNods(i,j,1,:),[1,2])*distNeighbNods(i,j,1);
    vec_P_N=reshape(uVecNeighbNods(i,j,2,:),[1,2])*distNeighbNods(i,j,2);
    
    %face west********************************************************
    % mass flux fW
    
    %if possitive node  P is upwind and node W is downwind//////
    %if negative  node  W is upwind and node P is downwind////////
    
    if fW>0 %POSITIVE MASS FLUX
        %node  P upwind, node W downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_W;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i,j-1);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i,j-1);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_w=0;   
        else
            r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_w_u=VanLeerGradientLimiter(r_w_u);
        
            %compute deferred correction contribution for source
            % in face W
            suDC_w=-0.5*fW*psi_w_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_w=0;
            
        else
            r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_w_v=VanLeerGradientLimiter(r_w_v);
        
            %compute deferred correction contribution for source
            % in face W
            svDC_w=-0.5*fW*psi_w_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  W upwind, node P downwind
        grad_phi_up_u=grad_u_W;
        grad_phi_up_v=grad_v_W;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_W;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j-1);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j-1);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_w=0;
            
        else
            r_w_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_w_u=VanLeerGradientLimiter(r_w_u);
            %compute deferred correction contribution for source
            % in face W
            suDC_w=-0.5*fW*psi_w_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_w=0;
            
        else
            r_w_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_w_v=VanLeerGradientLimiter(r_w_v);
            %compute deferred correction contribution for source
            % in face W
            svDC_w=-0.5*fW*psi_w_v*denom_v;
        end
    
    end
    
    %face north*******************************************************
    % mass flux fN
    %if possitive node  P is upwind and node N is downwind
    %if negative  node  N is upwind and node P is downwind
    
    if fN>0 %POSITIVE MASS FLUX
        %node  P upwind, node N downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_N;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i-1,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i-1,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_n=0;   
        else
            r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_n_u=VanLeerGradientLimiter(r_n_u);
        
            %compute deferred correction contribution for source
            % in face N
            suDC_n=-0.5*fN*psi_n_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_n=0;
            
        else
            r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_n_v=VanLeerGradientLimiter(r_n_v);
        
            %compute deferred correction contribution for source
            % in face N
            svDC_n=-0.5*fN*psi_n_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  N upwind, node P downwind
        grad_phi_up_u=grad_u_N;
        grad_phi_up_v=grad_v_N;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_N;
        %upwind and downwind phi values for u
        phi_up_u=u(i-1,j);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i-1,j);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_n=0;
            
        else
            r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_n_u=VanLeerGradientLimiter(r_n_u);
            %compute deferred correction contribution for source
            % in face N
            suDC_n=-0.5*fN*psi_n_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_n=0;
            
        else
            r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_n_v=VanLeerGradientLimiter(r_n_v);
            %compute deferred correction contribution for source
            % in face N
            svDC_n=-0.5*fN*psi_n_v*denom_v;
        end
    
    end
    
    %total deferred correction
    %DC for U
    suDC_total= suDC_w + suDC_n ;
    %DC for V
    svDC_total= svDC_w + svDC_n ;

    
    %Sources (Pressure dE(i,j)rivative)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + suDC_total;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
    
    %1.3.4- South-West (Symetric at South and inlet at West) *************
    i=ny;
    j=1;

    %get face areas
    lgt_fw=reshape(lgtFaces(i,j,1),[1 1]);
    lgt_fn=reshape(lgtFaces(i,j,2),[1 1]);
    lgt_fe=reshape(lgtFaces(i,j,3),[1 1]); 
    %lgt_fs=reshape(lgtFaces(i,j,4),[1 1]);
    
    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_face(i,j);
    fE=rho*lgt_fe*u_face(i,j+1); 
    fN=rho*lgt_fn*v_face(i,j); 
    %fS=-rho*lgt_fs*v_face(i+1,j);
    
    %coeffitients using upwind squeme
    aW(i,j)=dW(i,j);
    aE(i,j)=dE(i,j) + max(0,-fE);
    aN(i,j)=dN(i,j) + max(0,-fN);
    aS(i,j)=dS(i,j);%No applicable for u due to SBC
    aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j)  +fN +fE;
    aPv(j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j)  +fN +fE;

    %____________ CROSS DIFFUSION TERMS______________________

    %______________________U__________________________
            
    %Get velocity vertexes values
    vert_u=reshape(u_corners(i,j,:,:),[1,4]);

    %Split velocities
    u_wn=vert_u(1);
    u_en=vert_u(2);
    u_es=vert_u(3);
    u_ws=vert_u(4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values
    vert_v=reshape(v_corners(i,j,:,:),[1,4]);

    %Split velocities
    v_wn=vert_v(1);
    v_en=vert_v(2);
    v_es=vert_v(3);
    v_ws=vert_v(4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell
    gradPcell=reshape(grad_p(i,j,:,:),[1,2]);
    
    %Grad in x
    grad_px=gradPcell(1);
    %Grad in y
    grad_py=gradPcell(2);

    %____ TVD TOTAL CONTRIBUTION FOR THE CELL   ____________________

    %____________  TRANSPORTED VARIABLE = U   ______________________
            
    %Get gradients at central node and neighborhod nodes 
    
    grad_u_P=reshape(grad_u(i,j,:,:),[1,2]);
    grad_u_N=reshape(grad_u(i-1,j,:,:),[1,2]);
    grad_u_E=reshape(grad_u(i,j+1,:,:),[1,2]);
    
     %____________  TRANSPORTED VARIABLE = V    _____________________
                
    grad_v_P=reshape(grad_v(i,j,:,:),[1,2]);
    grad_v_N=reshape(grad_v(i-1,j,:,:),[1,2]);
    grad_v_E=reshape(grad_v(i,j+1,:,:),[1,2]);
    
    %vectors from central node_p to neighboorhod nodes (unitare vector
    %*magnitude)
    vec_P_N=reshape(uVecNeighbNods(i,j,2,:),[1,2])*distNeighbNods(i,j,2);
    vec_P_E=reshape(uVecNeighbNods(i,j,3,:),[1,2])*distNeighbNods(i,j,3);
    
    
    
    %face north*******************************************************
    % mass flux fN
    %if possitive node  P is upwind and node N is downwind
    %if negative  node  N is upwind and node P is downwind
    
    if fN>0 %POSITIVE MASS FLUX
        %node  P upwind, node N downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_N;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i-1,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i-1,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_n=0;   
        else
            r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_n_u=VanLeerGradientLimiter(r_n_u);
        
            %compute deferred correction contribution for source
            % in face N
            suDC_n=-0.5*fN*psi_n_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_n=0;
            
        else
            r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_n_v=VanLeerGradientLimiter(r_n_v);
        
            %compute deferred correction contribution for source
            % in face N
            svDC_n=-0.5*fN*psi_n_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  N upwind, node P downwind
        grad_phi_up_u=grad_u_N;
        grad_phi_up_v=grad_v_N;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_N;
        %upwind and downwind phi values for u
        phi_up_u=u(i-1,j);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i-1,j);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_n=0;
            
        else
            r_n_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_n_u=VanLeerGradientLimiter(r_n_u);
            %compute deferred correction contribution for source
            % in face N
            suDC_n=-0.5*fN*psi_n_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_n=0;
            
        else
            r_n_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_n_v=VanLeerGradientLimiter(r_n_v);
            %compute deferred correction contribution for source
            % in face N
            svDC_n=-0.5*fN*psi_n_v*denom_v;
        end
    
    end
    
    %face east********************************************************
    % mass flux fE
    %if possitive node P is upwind and node E is downwind
    %if negative  node E is upwind and node P is downwind
    
    if fE>0 %POSITIVE MASS FLUX
        %node  P upwind, node E downwind
        grad_phi_up_u=grad_u_P;
        grad_phi_up_v=grad_v_P;
        %vector from upwind node to downwind node 
        up_do_vec=vec_P_E;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j);
        phi_do_u=u(i,j+1);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j);
        phi_do_v=v(i,j+1);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
        
        %DEFERRED CORRECTION FOR U
        if denom_u==0 %DENOMINATOR IS ZERO
            suDC_e=0;   
        else
            r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_e_u=VanLeerGradientLimiter(r_e_u);
        
            %compute deferred correction contribution for source
            % in face E
            suDC_e=-0.5*fE*psi_e_u*denom_u;
            
        end
        %DEFERRED CORRECTION FOR V
        if denom_v==0
            svDC_e=0;
            
        else
            r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_e_v=VanLeerGradientLimiter(r_e_v);
        
            %compute deferred correction contribution for source
            % in face E
            svDC_e=-0.5*fE*psi_e_v*denom_v;
        end
    
    else %NEGATIVE MASS FLUX
        %node  E upwind, node P downwind
        grad_phi_up_u=grad_u_E;
        grad_phi_up_v=grad_v_E;
        %vector from upwind node to downwind node 
        up_do_vec=-vec_P_E;
        %upwind and downwind phi values for u
        phi_up_u=u(i,j+1);
        phi_do_u=u(i,j);
        %upwind and downwind phi values for v
        phi_up_v=v(i,j+1);
        phi_do_v=v(i,j);
        %denominators
        denom_u=phi_do_u -phi_up_u;
        denom_v=phi_do_v -phi_up_v;
    
        %DEFERRED CORRECTION FOR U
        if denom_u==0%DENOMINATOR IS ZERO
            suDC_e=0;
            
        else
            r_e_u=(2*dot(grad_phi_up_u,up_do_vec)/denom_u) - 1;
            psi_e_u=VanLeerGradientLimiter(r_e_u);
            %compute deferred correction contribution for source
            % in face E
            suDC_e=-0.5*fE*psi_e_u*denom_u;
        end
        
        %DEFERRED CORRECTION FOR v
        if denom_v==0%DENOMINATOR IS ZERO
            svDC_e=0;
            
        else
            r_e_v=(2*dot(grad_phi_up_v,up_do_vec)/denom_v)-1;
            psi_e_v=VanLeerGradientLimiter(r_e_v);
            %compute deferred correction contribution for source
            % in face E
            svDC_e=-0.5*fE*psi_e_v*denom_v;
        end
    
    end
    
    %total deferred correction
    %DC for U
    suDC_total= suDC_n + suDC_e ;
    %DC for V
    svDC_total= svDC_n + svDC_e ;

    %Sources (Pressure dE(i,j)rivative)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(aW(i,j) - fW) + suDC_total;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + svDC_total;
        
end

