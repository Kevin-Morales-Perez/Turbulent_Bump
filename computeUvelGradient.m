function [grad_u] = computeUvelGradient(u,wlsqOperator,grad_u,u0,...
    nx_upstr,nx_dwnstr)
%Computation of gradients at the cell centers

%pressure 
size_fields= size(u);
ny=size_fields(1);
nx=size_fields(2);


%internal cells 
for i= 2:ny-1
    for j=2:nx-1

        %get central and neighborhood nodes;
        u_w=u(i,j-1);
        u_n=u(i-1,j);
        u_e=u(i,j+1);
        u_s=u(i+1,j);
        u_p=u(i,j); 

        difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

        %Compute gradient
        grad_u(i,j,:,:)=(wlsqOp*difu_vec)';
       
    end
end

% IN EDGES AND CORNERS CONSIDER BOUNDARY CONDITIONS FOR U velocity!  _____

%________________________EDGES_________________________
%West Edge (
j=1;
for i=2:ny-1

    %get central and neighborhood nodes;

    u_w=u0; %Inlet
    u_n=u(i-1,j);
    u_e=u(i,j+1);
    u_s=u(i+1,j);
    u_p=u(i,j);

    difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_u(i,j,:,:)=(wlsqOp*difu_vec)';


end


%North Edge
i=1;
for j=2:nx-1

    %get central and neighborhood nodes;
    u_w=u(i,j-1);
    u_n=u(i,j); %Neumman condition
    u_e=u(i,j+1);
    u_s=u(i+1,j);
    u_p=u(i,j); 

    difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_u(i,j,:,:)=(wlsqOp*difu_vec)';
    
end


%East Edge
j=nx;
for i=2:ny-1
    
    %get central and neighborhood nodes;
    u_w=u(i,j-1);
    u_n=u(i-1,j);
    u_e=u(i,j);%outlet condition (Neumman)
    u_s=u(i+1,j);
    u_p=u(i,j); 

    difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_u(i,j,:,:)=(wlsqOp*difu_vec)';


end


%South Edge
i=ny;
for j=2:nx-1

    if (j<=nx_upstr) || j>(nx-nx_dwnstr)

        %get central and neighborhood nodes;
        u_w=u(i,j-1);
        u_n=u(i-1,j);
        u_e=u(i,j+1);
        u_s=u(i,j);% Symetric condition
        u_p=u(i,j); 
        
        difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec
        
        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);
        
        %Compute gradient
        grad_u(i,j,:,:)=(wlsqOp*difu_vec)';

    else
        %over the plate
        %get central and neighborhood nodes;
        u_w=u(i,j-1);
        u_n=u(i-1,j);
        u_e=u(i,j+1);
        u_s=0;% No slip condition
        u_p=u(i,j); 
        
        difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec
        
        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);
        
        %Compute gradient
        grad_u(i,j,:,:)=(wlsqOp*difu_vec)';

    end


end


%______________________Corners______________________________
%west north
i=1;
j=1;

%get central and neighborhood nodes;
u_w=u0; %Inlet
u_n=u(i,j); %Neumman condition
u_e=u(i,j+1);
u_s=u(i+1,j);
u_p=u(i,j); 

difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_u(i,j,:,:)=(wlsqOp*difu_vec)';


%east north
i=1;
j=nx;

%get central and neighborhood nodes;

u_w=u(i,j-1);
u_n=u(i,j); %Neumman condition
u_e=u(i,j);%outlet condition
u_s=u(i+1,j);
u_p=u(i,j);

difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_u(i,j,:,:)=(wlsqOp*difu_vec)';


%east soth
i=ny;
j=nx;

%get central and neighborhood nodes;
u_w=u(i,j-1);
u_n=u(i-1,j);
u_e=u(i,j);%outlet condition
u_s=u(i,j);% Symetric condition
u_p=u(i,j); 

difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_u(i,j,:,:)=(wlsqOp*difu_vec)';


%west south
i=ny;
j=1;
 
%get central and neighborhood nodes;
u_w=u0; %Inlet
u_n=u(i-1,j);
u_e=u(i,j+1);
u_s=u(i,j);% Symetric condition
u_p=u(i,j); 

difu_vec=[u_w;u_n;u_e;u_s]-u_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_u(i,j,:,:)=(wlsqOp*difu_vec)';
