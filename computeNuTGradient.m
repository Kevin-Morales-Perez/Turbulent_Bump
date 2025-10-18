function [grad_nu_tilde] = computeNuTGradient(nu_tilde,wlsqOperator,...
    grad_nu_tilde,nu,nx_upstr,nx_dwnstr)
%Computation of gradients at the cell centers

%nu_tilde 
size_fields= size(nu_tilde);
ny=size_fields(1);
nx=size_fields(2);


%internal cells 
for i= 2:ny-1
    for j=2:nx-1

        %get central and neighborhood nodes;
        nut_w=nu_tilde(i,j-1);
        nut_n=nu_tilde(i-1,j);
        nut_e=nu_tilde(i,j+1);
        nut_s=nu_tilde(i+1,j);
        nut_p=nu_tilde(i,j); 

        difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

        %Compute gradient
        grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';
       
    end
end

% IN EDGES AND CORNERS CONSIDER BOUNDARY CONDITIONS FOR NU TILDE!  _____

%________________________EDGES_________________________
%West Edge (
j=1;
for i=2:ny-1

    %get central and neighborhood nodes;

    nut_w=0.1*nu; %Inlet
    nut_n=nu_tilde(i-1,j);
    nut_e=nu_tilde(i,j+1);
    nut_s=nu_tilde(i+1,j);
    nut_p=nu_tilde(i,j);

    difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';


end


%North Edge
i=1;
for j=2:nx-1

    %get central and neighborhood nodes;
    nut_w=nu_tilde(i,j-1);
    nut_n=nu_tilde(i,j); %Neumman condition
    nut_e=nu_tilde(i,j+1);
    nut_s=nu_tilde(i+1,j);
    nut_p=nu_tilde(i,j); 

    difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';
    
end


%East Edge
j=nx;
for i=2:ny-1
    
    %get central and neighborhood nodes;
    nut_w=nu_tilde(i,j-1);
    nut_n=nu_tilde(i-1,j);
    nut_e=nu_tilde(i,j);%outlet condition (Neumman)
    nut_s=nu_tilde(i+1,j);
    nut_p=nu_tilde(i,j); 

    difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';


end


%South Edge
i=ny;
for j=2:nx-1

    if (j<=nx_upstr) || j>(nx-nx_dwnstr)

        %get central and neighborhood nodes;
        nut_w=nu_tilde(i,j-1);
        nut_n=nu_tilde(i-1,j);
        nut_e=nu_tilde(i,j+1);
        nut_s=nu_tilde(i,j);% Symetric condition
        nut_p=nu_tilde(i,j); 
        
        difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec
        
        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);
        
        %Compute gradient
        grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';

    else
        %over the plate
        %get central and neighborhood nodes;
        nut_w=nu_tilde(i,j-1);
        nut_n=nu_tilde(i-1,j);
        nut_e=nu_tilde(i,j+1);
        nut_s=0;% No eddy viscosity at solid surface
        nut_p=nu_tilde(i,j); 
        
        difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec
        
        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);
        
        %Compute gradient
        grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';

    end


end


%______________________Corners______________________________
%west north
i=1;
j=1;

%get central and neighborhood nodes;
nut_w=0.1*nu; %Inlet
nut_n=nu_tilde(i,j); %Neumman condition
nut_e=nu_tilde(i,j+1);
nut_s=nu_tilde(i+1,j);
nut_p=nu_tilde(i,j); 

difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';


%east north
i=1;
j=nx;

%get central and neighborhood nodes;

nut_w=nu_tilde(i,j-1);
nut_n=nu_tilde(i,j); %Neumman condition
nut_e=nu_tilde(i,j);%outlet condition
nut_s=nu_tilde(i+1,j);
nut_p=nu_tilde(i,j);

difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';


%east soth
i=ny;
j=nx;

%get central and neighborhood nodes;
nut_w=nu_tilde(i,j-1);
nut_n=nu_tilde(i-1,j);
nut_e=nu_tilde(i,j);%outlet condition
nut_s=nu_tilde(i,j);% Symetric condition
nut_p=nu_tilde(i,j); 

difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';


%west south
i=ny;
j=1;
 
%get central and neighborhood nodes;
nut_w=0.1*nu; %Inlet
nut_n=nu_tilde(i-1,j);
nut_e=nu_tilde(i,j+1);
nut_s=nu_tilde(i,j);% Symetric condition
nut_p=nu_tilde(i,j); 

difnut_vec=[nut_w;nut_n;nut_e;nut_s]-nut_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_nu_tilde(i,j,:,:)=(wlsqOp*difnut_vec)';
