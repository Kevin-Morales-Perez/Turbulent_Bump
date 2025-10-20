%FVM BUMP IN CHANNEL WITH INFLATION LAYER IN BOUNDARY
%LEADING AND TRAILING EDGE
%FREE FLOW UPSTREAM AND DOWNSTREAM ZONES
%ALGORITHM: SIMPLE
%SCHEME: UPWIND
% MATRIX SOL. METHOD: GAUSS
%FLUID: AIR AT 20° AND SEA LEVEL
%FACE INTERPOLATION METHOD: RIE-CHOW
%CO-LOCATED NON-ORTHOGONAL GRID 
%STEADY STATE INCOMPRESSIBLE REYNOLDS AVERAGED NAVIER STOKES
%TURBULENCE MODEL: SPALART - ALLMARAS 

%BOUNDARY  CONDITIONS

%        BOUNDARY CONDITIONS FOR FLOW OVER A BUMP
%              NEUMMAN  (DU/DY=0 ,DV/DY=0)
%                                                      O N  
%I                                                     U E   
%N                                                     T U  DU/DX=0
%L                                                     L M
%E                                                     E M  DV/DX=0
%T                                                     T A 
% ___ SYMETRIC ___|___   NO SLIP ___|___ SYMETRIC  ____  N



close all

%tic


%% DATA

%Geometry_________________________________________________________________
bumpLgt=1.5;%Match with MeshBump.m

%Fluid____________________________________________________________________
rho=1.2;                    %Density (Kg/m3)
mu=1.8e-5;                  %Dynamic Viscosity (N*s/m^2)
nu = mu/rho;                 %Kinematic Viscosity (m/s^2)
u0=0.0125;                   %Velocity at the inlet (m/s) 
p0=1;                       %Outlet pressure (Prescribed)
Re=u0*rho*(bumpLgt-0.5)/mu; %Reynolds number

% 0.0771; 45

%Spalart Allmaras Model Constants_________________________________________
%General constans
kappa=0.41; %Karman Constant
sigma_sa=2/3;%Turbulent Prantl Number used in SA model
%Basic constants
cb1=0.1355;
cb2=0.622;
%Viscous Damping
cv1=7.1;
%Wall Turbulent Destruction 
cw1=cb1/(kappa^2) +  (1 + cb2)/sigma_sa;
cw2=0.3;
cw3=2;

%% GRID 

%run MeshBump.m  to generate the mesh

%Geometry defined in MeshBump.m

run("meshbump2.m") %<-----------------------------------------------

%Create variables to store geometrical properties 

Xctrs=zeros(ny,nx);%Mesh with cell center positions in X axis

Yctrs=zeros(ny,nx);%Mesh with cell center positions in Y axis 

cellVertxs=zeros(ny,nx,4,2);%Array with coordinates of vertex of each
%  cells [Vwn;Ven;Ves;Vws]=reshape(cellVertxs(i,j,:,:),[4,2])*

cellCentrs=zeros(ny,nx,2);%Matrix with centroids of each cell;
%[xi,ji]=reshape(cellCentrs(i,j,:),[1,2])*

lgtFaces=zeros(ny,nx,4); %matrix which contains lenght of the 4 faces of 
% each cell ,[lgtFw,lgtFn,lgtFe,lgtFs]=reshape(lgtFaces(i,j,:),[1 4])*

faceCentrs=zeros(ny,nx,4,2);%Array which contains the coordinates of
% points at the centres of each face
%[fCw;fCn,fCe,fCs]=reshape(faceCentrs(i,j,:,:),[4,2])*

cellVols=zeros(ny,nx);%matrix which contains the volume of each cell

uVecNormFaces=zeros(ny,nx,4,2);% Array with unitary vectors normal
%to the faces ,[nw;nn;ne;ns]=reshape(uVecNormFaces(i,j,:,:),[4,2])
%In Malalasekera it is defined as n*

uVecNeighbNods=zeros(ny,nx,4,2); %Array with unitary vectors from central
%node to Neighborhood nodes 
% [uNnw;uNnn;uNne;uNns]=reshape(uVecNeighbNods(i,j,:,:),[4,2])
%In Malalasekera it is defined as Xi
%if there is no neighborhood node and instead is a boundary the node is
%conected to the face center 

distNeighbNods=zeros(ny,nx,4);%Array with the distance between central 
%node and neigborhood nodes (If there is no neighborhood node is calculated
%the distance to the center of the face instead
%[distNw,distNw,distNw,distNw]=reshape(distNeighbNods(i,j,:,:),[1,4])

uVecParlFaces=zeros(ny,nx,4,2);%Array wich contains unitary vectors 
%paralel to each face of the cell (All from weast to east and from s to n)
%[nw;nn;ne;ns]=reshape(uVecParlFaces(i,j,:,:),[4,2])
%In malalasekera it is defined as eta *

wlsqOperator=zeros(ny,nx,2,4);%Array for weighted least squares 
% operators 
%wlsop=reshape(wlsqOperator(i,j,:,:),[2,4])
%(wlsop*dif_vec)' ([2,4][4,1])'=[1,2]

VecCentNodFaceCents=zeros(ny,nx,4,2);%Array that contains vectors from 
%cell center to each face center 
%[vpw;vpn;vpe;vps]=reshape(VecCentNodFaceCents(i,j,:,:),[4,2])

VecCentVertx=zeros(ny,nx,4,2);% Array that contains vectors from cell 
%center to each vertex of the cell 
%[vc_wn;vc_en;vc_es;vc_ws;]=reshape(VecCentVertx(i,j,:,:),[4,2])

weightDistFactors=zeros(ny,nx,4);% Weight distance factor to
%calculate influence of cell center in face center
%reshape(weightDistFactors(i,j,:),[1,4])

distMinWall=zeros(ny,nx);%Minimun distance to the nearest wall for Spalart -
%  Allmaras Equation (Will be used further)

%process the mesh to calculate geometrical properties
[cellVertxs,cellCentrs,lgtFaces,uVecParlFaces,faceCentrs,...
    cellVols,uVecNormFaces,distNeighbNods,uVecNeighbNods,...
    wlsqOperator,VecCentNodFaceCents,VecCentVertx,weightDistFactors,...
    Xctrs,Yctrs,distMinWall] =mesh_geometrical_process(X,Y,cellVertxs,...
    cellCentrs,lgtFaces,uVecParlFaces,faceCentrs,cellVols,uVecNormFaces,...
    distNeighbNods,uVecNeighbNods,wlsqOperator,VecCentNodFaceCents,...
    VecCentVertx,weightDistFactors,Xctrs,Yctrs,distMinWall,distPlat);

%% VARIABLES AND COEFFITIENS

u=zeros(ny,nx); %Velocity in X axis
v=zeros(ny,nx);% Velocity in Y axis 
p=ones(ny,nx);%Pressure

% Reynolds stress tensor  Tau_ij
tau_xx=zeros(ny,nx); %-rho*u'2- Normal
tau_xy=zeros(ny,nx); %-rho*u'v'- Shear
tau_yy=zeros(ny,nx); %-rho*v'2- Normal

% Eddy viscosity
mu_turbulent=zeros(ny,nx); %Eddy viscosity
nu_tilde=zeros(ny,nx) + 0.1*nu; % Transported viscosity in Spalart - 
% Allmaras model

%velocities normal to faces
u_face=zeros(ny,nx+1);%for faces W and E
v_face=zeros(ny+1,nx);%for faces N and S

%velocities at corners for cross - gradient computation
u_corners=zeros(ny,nx,1,4);%u
v_corners=zeros(ny,nx,1,4);%v

%Vorticity
vorticity=zeros(ny,nx);

%Gradient of velocity and pressure 
grad_u=zeros(ny,nx,1,2);
grad_v=zeros(ny,nx,1,2);
grad_p=zeros(ny,nx,1,2);
grad_p_prime=zeros(ny,nx,1,2);

%Gradient of Reynolds Stresses
grad_tau_xx=zeros(ny,nx,1,2);
grad_tau_xy=zeros(ny,nx,1,2);
grad_tau_yy=zeros(ny,nx,1,2);

%Grad of nu_tilde 
grad_nu_tilde=zeros(ny,nx,1,2);

%Nu~ at corners
nu_tilde_corners=zeros(ny,nx,1,4);

%Grad of Eddy Viscosity
grad_mu_turbulent=zeros(ny,nx,1,2);
%Eddy viscosity at faces
mu_turbulent_fw=zeros(ny,nx);%Eddy viscosity at face W
mu_turbulent_fn=zeros(ny,nx);%Eddy viscosity at face N
mu_turbulent_fe=zeros(ny,nx);%Eddy viscosity at face E
mu_turbulent_fs=zeros(ny,nx);%Eddy viscosity at face S

%reshape(grad_u(i,j,:,:),[1,2])
u_face(:,1)=u0;%Imposing initial velocity at the inlet 

%Guessed nodal pressure and velocity fields
u_star=zeros(ny,nx); %X component
v_star=zeros(ny,nx); %Y component

%Correction pressure  field
p_prime = zeros(ny,nx); %Pressure 

%Momentum equation coeffitients
aW=zeros(ny,nx); %Weast
aN=zeros(ny,nx); %North
aE=zeros(ny,nx); %East
aS=zeros(ny,nx); %South
aP=zeros(ny,nx); %Diagonal
aPv=zeros(1,nx);%Diagonal Coeffitienst for simetric boundary condition
%at the free-stream zones

%Pressure equation coeffitients
ap_W=zeros(ny,nx); %Weast
ap_N=zeros(ny,nx); %North
ap_E=zeros(ny,nx); %East
ap_S=zeros(ny,nx); %South
ap_P=zeros(ny,nx); %Diagonal

%Nu_tilde equation coeffitients
ant_W=zeros(ny,nx); %Weast
ant_N=zeros(ny,nx); %North
ant_E=zeros(ny,nx); %East
ant_S=zeros(ny,nx); %South
ant_P=zeros(ny,nx); %Diagonal

%Source terms 
suX=zeros(ny,nx); %Source for U
suY=zeros(ny,nx); %Source for V
suP=zeros(ny,nx); %Source for P
suNt=zeros(ny,nx); %Source for Nu_tilde 

%Residual
rsid_x=zeros(ny,nx); %X momentum
rsid_y=zeros(ny,nx); %Y momentum
rsid_p=zeros(ny,nx);% Pressure correction eq.
rsid_cont=zeros(ny,nx); % Continuity
rsid_nt=zeros(ny,nx); %nu_tilde 


%error from residual
err_x=1;
err_y=1;
err_p=1;
err_nt=1;
%main iterations counter 
iterations_cont=0;

%Error for interior iterations
epsilon_uv=0.9e-20;
epsilon_p=5e-9;
epsilon_nt=9e-9;

%Underelaxation factors
alpha_uv=0.3;  %<---------------------     x-y momentum 
alpha_p=0.0008; %<---------------------     pressure correction
alpha_uv2=1;    %<---------------------     Rie - chow face velocity


max_iterations=10000;% Max outer iterations <---------------------
max_iterations_u=120;% Max iterations for momentum eqs.
max_iterations_v=120;% Max iterations for momentum  y eq .
max_iterations_p=120;% Max iterations for pressure eq.
max_iterations_nt=1;%Max iterations for nu_tilde SA Transport eq. 
residualsMat=zeros(max_iterations,4);% Residual Matrix
error_tgt=1.15e-18; %Target error 5.1e-18 for Re=2000
max_residual=1e10; % unstable value flag
convergedFlg=false; %Flag for convergence 


%%  COMPUTATION OF DIRECT AND NON ORTHOGONAL DIFFUSSION COEFFITIENTS
%Diffusive fluxes (Direct Gradient terms)
dW=zeros(ny,nx); %Face W
dN=zeros(ny,nx); %Face N
dE=zeros(ny,nx); %Face E
dS=zeros(ny,nx); %Face S

%Cross difussion terms (but they will be added as a source)
dW_c=zeros(ny,nx); %Face W
dN_c=zeros(ny,nx); %Face N
dE_c=zeros(ny,nx); %Face E
dS_c=zeros(ny,nx); %Face S

%compute coeffitients
[dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c] = ...
    diffusive_coeffitients(nx,ny,distNeighbNods,uVecNormFaces,...
    uVecNeighbNods,lgtFaces,uVecParlFaces,1,dW,dN,dE,dS,dW_c,dN_c,...
    dE_c,dS_c);

%% SOLVER 

tic

% Main loop
while convergedFlg==false
    %Iterations counter
    iterations_cont=iterations_cont+1;
    
    %Reset continuity residual for the new iteration
    rsid_cont(:) = 0;

        %0.- COMPUTE PRESSURE  GRADIENT

    grad_p = computePressGradient(p,wlsqOperator,grad_p,p0);
    
    %Compute velocity at vertex at each new step for non orthogonal
    %Difussion terms (Can be Simplified)

    u_corners = computePhiVertex(u,grad_u,VecCentVertx,u_corners);
    v_corners= computePhiVertex(v,grad_v,VecCentVertx,v_corners);

    %compute eddy viscosity Gradient
    [grad_mu_turbulent] = computeMuTurbGradient(mu_turbulent,...
    wlsqOperator,grad_mu_turbulent,mu,nx_upstr,nx_dwnstr);

    %compute eddy viscosity at faces;

    [mu_turbulent_fw,mu_turbulent_fn,mu_turbulent_fe,...
    mu_turbulent_fs] = computePhiFaces(grad_mu_turbulent,mu_turbulent,...
    VecCentNodFaceCents,nx,ny,mu_turbulent_fw,mu_turbulent_fn,...
    mu_turbulent_fe,mu_turbulent_fs);


    %_____________________SIMPLE ALGORITHM______________________________
    
    %1.- MOMENTUM LINK COEFFITIENTS AND SOURCES
    %[aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeff(nx,ny,...
    %nx_upstr,nx_dwnstr,rho,lgtFaces,u_face,v_face,u0,...
    %aW,aE,aN,aS,dW,dE,dN,dS,dW_c,dN_c,dE_c,dS_c,aP,aPv,suX,suY,...
    %cellVols,u_corners,v_corners,grad_p,mu,mu_turbulent_fw,...
    %mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs);

    [aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeffTVD(nx,ny,...
    nx_upstr,nx_dwnstr,rho,lgtFaces,u_face,v_face,u0,...
    aW,aE,aN,aS,dW,dE,dN,dS,dW_c,dN_c,dE_c,dS_c,aP,aPv,suX,suY,...
    cellVols,u_corners,v_corners,grad_p,mu,mu_turbulent_fw,...
    mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs,grad_u,grad_v,...
    uVecNeighbNods,distNeighbNods,u,v);

    %2.- SOLVE X MOMENTUM

    [u,rsid_x,err_x] = solve_momentum(err_x,max_iterations_u,...
    nx,ny,u,alpha_uv,aP,aW,aN,aE,aS,suX,u_star,rsid_x);

    %3.- SOLVE Y MOMENTUM

    aP_temp_v=aP;% temporal central coeffitients matrix for v

    aP_temp_v(ny,nxSymcond)=aPv(nxSymcond);%Replace coeffitients that 
    % should not be the same for the Simetryc boundary condition

    [v,rsid_y,err_y] = solve_momentum(err_y,max_iterations_v,...
    nx,ny,v,alpha_uv,aP_temp_v,aW,aN,aE,aS,suY,v_star,rsid_y);

    aP=aP/alpha_uv;
    aPv=aPv/alpha_uv;

    %4.- FACE VELOCITY COMPUTATION USING RIE - CHOW INTERPOLATION*
 
    [u_face,v_face] = face_vel_intRCnonrthgn2(u,v,p,u_face,v_face,...
    uVecNormFaces,distNeighbNods,uVecNeighbNods,weightDistFactors,...
    grad_p,aP,aPv,cellVols,nx_upstr,nx_dwnstr);

    %Neumman boundary condition for outlet flow at east edge and north edge
    u_face(:,end)=u(:,end);
    v_face(1,:)=v(1,:);

    %5.- PRESSURE CORRECTION LINK COEFFITIENTS AND MASS INBALANCE (SOURCE)

    [ap_W,ap_N,ap_E,ap_S,ap_P,suP] = pressure_link_coeff(u_face,...
    v_face,aP,aPv,ap_W,ap_N,ap_E,ap_S,ap_P,suP,weightDistFactors,...
    cellVols,distNeighbNods,lgtFaces,nx_upstr,nx_dwnstr);

    %6.- SOLVE PRESSURE CORRECTION

    % Reset p_prime for the new iteration
    p_prime=zeros(ny,nx);

    [p_prime,rsid_p,err_p] = solve_presscorr(err_p,max_iterations_p...
        ,nx,ny,p_prime,ap_P,ap_W,ap_N,ap_E,ap_S,suP,rsid_p...
        ,epsilon_p);

    %7.- CORRECT PRESSURE

    p=p + alpha_p*p_prime;

    %8.- CORRECT CELL CENTER VELOCITY

        %8.0 - COMPUTE PRESSURE CORRECTION GRADIENT
    [grad_p_prime] = computepPrimeGradient(p_prime,wlsqOperator,...
        grad_p_prime);

    [u_star,v_star] = cvel_correct(aP,aPv,u,v,grad_p_prime,cellVols,...
    nx_upstr,nx_dwnstr,alpha_uv2);

    %9.- CORRECT FACE VELOCITY
    
    [u_face,v_face] = fvel_correct(u_face,v_face,p_prime,aP,aPv,...
    cellVols,weightDistFactors,distNeighbNods,nx_upstr,nx_dwnstr,alpha_uv2);

    %______________________TURBULENCE MODELLING__________________________

   
    %compute velocity gradients
    grad_u= computeUvelGradient(u,wlsqOperator,grad_u,u0,nx_upstr,...
       nx_dwnstr);
    grad_v = computeVvelGradient(v,wlsqOperator,grad_v);
    %compute vorticity
    [vorticity] = computationVorticity(grad_u,grad_v,vorticity,ny,nx);
     
    %compute nu~ gradient
    grad_nu_tilde=computeNuTGradient(nu_tilde,wlsqOperator,...
        grad_nu_tilde,nu,nx_upstr,nx_dwnstr);
    
    %compute nu~ at vertexes

    nu_tilde_corners= computePhiVertex(nu_tilde,grad_nu_tilde,...
        VecCentVertx,nu_tilde_corners);

    %10.- SPALART -ALLMARAS MODEL LINK COEFFITIENTS AND SOURCE

    [ant_W,ant_N,ant_E,ant_S,ant_P,suNt] = saTurbulence_link_coeff(...
    nu_tilde,nu,vorticity,lgtFaces,cellVols,u_face,v_face,distMinWall,...
    grad_nu_tilde,dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c,nu_tilde_corners,...
    ant_W,ant_N,ant_E,ant_S,ant_P,suNt,...
    kappa,sigma_sa,cb1,cb2,cv1,cw1,cw2,cw3,...
    nx,ny,nxSolid);

    %11.-  SOLVE NU- TILDE

    [nu_tilde,rsid_nt,err_nt] = solve_presscorr(err_nt,max_iterations_nt...
        ,nx,ny,nu_tilde,ant_P,ant_W,ant_N,ant_E,ant_S,suNt,rsid_nt...
        ,epsilon_p);

    %[nu_tilde,rsid_nt,err_nt]=Gauss_general_solver(nu_tilde,nu_tilde,...
    %    ant_P,ant_W,ant_N,ant_E,ant_S,suNt,max_iterations_nt,rsid_nt,...
    %    epsilon_nt,alpha_sa);

    %12.- COMPUTE EDDY VISCOSITY FROM NU TILDE 

    [mu_turbulent] = saEddyViscosity(nu_tilde,nu,rho,cv1);
    
   
    %13.- CHECK RESIDUALS

    %fprintf('%.18f',residualsMat(iterations_cont,3))
    %fprintf('%.15f',residualsMat(iterations_cont,4))
    %residualsMat(:,3)=abs(residualsMat(:,3)-0.000001776715282356);
    %residualsMat(:,4)=abs(residualsMat(:,4)-0.000001776715282356);
    continuity_scalF=0;
    err_cont=abs(rms(suP(:)) - continuity_scalF);

    residualsMat(iterations_cont,:)=[err_x,err_y,err_cont,err_nt];
    disp(residualsMat(iterations_cont,:)) %display current error
    
    if err_x < error_tgt && err_y < error_tgt && err_cont < error_tgt ...
            && err_nt < error_tgt
       convergedFlg=true;
       fprintf("Converged  at iteration\n")
       disp(iterations_cont)
       save('flatplateFields.mat', 'u', 'v' , 'p',"nu_t" );
    elseif err_x>max_residual || err_y > max_residual ||...
            err_cont > max_residual || err_nt>max_residual
        fprintf("Unstable solution iterations stopped at iteration ")
        disp(iterations_cont)
        break
    elseif iterations_cont >= max_iterations
        fprintf("Max iterations reached \n")
        disp(iterations_cont)
        break
    elseif isnan(err_x) || isnan(err_y) || isnan(err_cont) || isnan(err_nt)
        fprintf("The system became undetermined  at iteration \n")
        disp(iterations_cont)
        break
    else
        convergedFlg=false;
    end
end


toc
beep()
convergedFlg=true;

%% POST PROCESS

fprintf('%s%d\n','Reynolds number ', Re)
fprintf('Total cells: %d\n', ncells);


%RESIDUALS 
figure(2)
plot(1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,1),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,2),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,3),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,4))
legend("U velocity","V velocity","Continuity"," Nu tilde")
title("Residuals")
xlabel("Iterations")
ylabel("Residual")
yscale log


if convergedFlg==true
    
    %COMPUTATION OF OTHER VARIABLES OF INTEREST

    %VELOCITY MAGNITUDE
    velocity_magnitude = sqrt(u.^2 + v.^2); 

    %PRESSURE COEFFITIENT
    x_pressure_coeff=zeros(size(nxSolid));

    for j=1:nx_fp
        j_1=nxSolid(j);
        x_pressure_coeff(j)=reshape(faceCentrs(ny,j_1,4,1),[1,1]);
    end 
    x_pressure_coeff=x_pressure_coeff-(domLgt-bumpLgt)/2;
    pressure_coeff=(p(ny,nxSolid)-p0)/rho*u0;

    %Compute gradient of Reynolds stresses
    grad_tau_xx=computeTauGradient(tau_xx,wlsqOperator,grad_tau_xx);
    grad_tau_xy=computeTauGradient(tau_xy,wlsqOperator,grad_tau_xy);
    grad_tau_yy=computeTauGradient(tau_yy,wlsqOperator,grad_tau_yy);

    %13.- BOUSSINESQ ASSUMPTION FOR SA MODEL 
    
    %compute turbulent stresses
    [tau_xx,tau_xy,tau_yy] = saBoussinesqTurbulentStresses(...
        mu_turbulent,grad_u,grad_v,tau_xx,tau_xy,tau_yy,nx,ny);
        

    figure(3)
    contourf(Xctrs,Yctrs,u, 20, 'LineColor', 'none')
    title("Velocity in x axis (m/s)")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal
    
    figure(4)
    contourf(Xctrs,Yctrs,v, 20, 'LineColor', 'none')
    title("Velocity in y axis (m/s)")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal
    
    
    figure(5)
    contourf(Xctrs, Yctrs, velocity_magnitude, 20, 'LineColor', 'none')
    colorbar
    %hold on 
    %quiver(Xctrs,Yctrs,u,v,'k')
    title("Velocity vector field with magnitude m/s")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    axis equal
    %hold off
    
    figure(6)
    contourf(Xctrs,Yctrs,p, 20, 'LineColor', 'none')
    title("pressure N/m^2")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal

    figure(7)
    streamslice(Xctrs, Yctrs, u, v); % For a sparse representation
    title("streamlines")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    axis equal

    
    figure(8)
    contourf(Xctrs, Yctrs, vorticity, 20, 'LineColor', 'none');
    colorbar
    title('Vorticity Field')
    xlabel('X')
    ylabel('Y')
    axis equal
    colormap jet;
    
    figure(9)
    u_prof=[0,flip(u(:,nx_upstr + 5)'),u0];
    ycord_prof=[0, flip(Yctrs(:,nx_upstr + 5)'),domHgt];
    plot(u_prof,ycord_prof,'- o')
    title("U profile of Flat Plate")
    xlabel("Y position")
    ylabel("U vel")

    figure(10)
    contourf(Xctrs,Yctrs,cellVols, 20, 'LineColor', 'none')
    title("Cell volumes")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal

    figure(11)
    plot(x_pressure_coeff,pressure_coeff ,'- o')
    title("PRESSURE COEFFITIENT ON THE PLATE ")
    xlabel("X POSITION")
    ylabel("PRESSURE COEFFITIENT ")

    
    figure(12)
    contourf(Xctrs,Yctrs,nu_tilde, 20, 'LineColor', 'none')
    title("Spalart - Allmaras nu~ transported variable")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal

    figure(13)
    contourf(Xctrs,Yctrs,mu_turbulent, 20, 'LineColor', 'none')
    title("Eddy Viscosity")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal

    figure(14)
    contourf(Xctrs,Yctrs,tau_xx, 20, 'LineColor', 'none')
    title("Turbulent Stress Tau_xx")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal

    figure(15)
    contourf(Xctrs,Yctrs,tau_xy, 20, 'LineColor', 'none')
    title("Turbulent Stress Tau_xy")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal

    figure(16)
    contourf(Xctrs,Yctrs,tau_yy, 20, 'LineColor', 'none')
    title("Turbulent Stress Tau_yy")
    xlabel("Lenght (m)")
    ylabel("Height (m)")
    colormap jet
    colorbar
    axis equal

end

figure(17)
contourf(Xctrs,Yctrs,suP, 20, 'LineColor', 'none')
title("Continuity Residual (Kg*m^2/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colorbar 
colormap jet
axis equal

%}







