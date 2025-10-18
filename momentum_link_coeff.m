function [aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeff(nx,ny,...
    nx_upstr,nx_dwnstr,rho,lgtFaces,u_face,v_face,u0,...
    aW,aE,aN,aS,dW,dE,dN,dS,dW_c,dN_c,dE_c,dS_c,aP,aPv,suX,suY,...
    cellVols,u_corners,v_corners,grad_p,mu,mu_turbulent_fw,...
    mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs)
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

            %_____ Sources (Pressure derivative + cross difussion) ______

            suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T;
            suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
    
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

        %Sources (Pressure derivative)
        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(aW(i,j) -fW);  
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
    
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


        %Sources (Pressure dE(i,j)rivative + boundary condition)
        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*aN(i,j); 
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
    
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


        %_____ Sources (Pressure derivative + cross difussion) ______

        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T;
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
    
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

          
            %_____ Sources (Pressure derivative + cross difussion) ______

            suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T;
            suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;

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

            suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T;
            suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;

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

    %Sources (Pressure dErivative + boundary condition)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(aW(i,j) + aN(i,j) -fW);
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
    
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

    %Sources (Pressure dE(i,j)rivative + boundary condition)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*aN(i,j);
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
    
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

    
    %Sources (Pressure dE(i,j)rivative)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
    
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

    %Sources (Pressure dE(i,j)rivative)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(aW(i,j) - fW);
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T;
        
end

