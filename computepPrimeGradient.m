function [grad_p_prime] = computepPrimeGradient(p_prime,wlsqOperator,...
    grad_p_prime)
%Computation of gradients at the cell centers

%pressure 
size_fields= size(p_prime);
ny=size_fields(1);
nx=size_fields(2);


%internal cells 
for i= 2:ny-1
    for j=2:nx-1

        %get central and neighborhood nodes;
        p_w=p_prime(i,j-1);
        p_n=p_prime(i-1,j);
        p_e=p_prime(i,j+1);
        p_s=p_prime(i+1,j);
        p_p=p_prime(i,j); 

        difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

        %Compute gradient
        grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';
       
    end
end

% IN EDGES AND CORNERS CONSIDER BOUNDARY CONDITIONS FOR PRESSURE CORRECTION
% !  _____

%________________________EDGES_________________________
%West Edge (
j=1;
for i=2:ny-1

    %get central and neighborhood nodes;

    p_w=0; %Inlet (velocity fixed , no pressure correction)
    p_n=p_prime(i-1,j);
    p_e=p_prime(i,j+1);
    p_s=p_prime(i+1,j);
    p_p=p_prime(i,j);

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';


end


%North Edge
i=1;
for j=2:nx-1

    %get central and neighborhood nodes;
    p_w=p_prime(i,j-1);
    p_n=p_prime(i,j);%Fixed velocity No pressure correction
    p_e=p_prime(i,j+1);
    p_s=p_prime(i+1,j);
    p_p=p_prime(i,j); 

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';
    
end


%East Edge
j=nx;
for i=2:ny-1
    
    %get central and neighborhood nodes;
    p_w=p_prime(i,j-1);
    p_n=p_prime(i-1,j);
    p_e=p_prime(i,j);%outlet condition, correction equal to cell center
    p_s=p_prime(i+1,j);
    p_p=p_prime(i,j); 

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';


end


%South Edge
i=ny;
for j=2:nx-1

    %get central and neighborhood nodes;
    p_w=p_prime(i,j-1);
    p_n=p_prime(i-1,j);
    p_e=p_prime(i,j+1);
    p_s=0;% Fixed velocity , no pressure correction needed
    p_p=p_prime(i,j); 

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';

end


%______________________Corners______________________________
%west north
i=1;
j=1;

%get central and neighborhood nodes;
p_w=0; %Inlet
p_n=p_prime(i,j);%Dirichlet condition
p_e=p_prime(i,j+1);
p_s=p_prime(i+1,j);
p_p=p_prime(i,j); 

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';


%east north
i=1;
j=nx;

%get central and neighborhood nodes;

p_w=p_prime(i,j-1);
p_n=p_prime(i,j);%Neummancondition
p_e=p_prime(i,j);%outlet condition
p_s=p_prime(i+1,j);
p_p=p_prime(i,j);

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';


%east soth
i=ny;
j=nx;

%get central and neighborhood nodes;
p_w=p_prime(i,j-1);
p_n=p_prime(i-1,j);
p_e=p_prime(i,j);%outlet condition
p_s=0;% No perpendicular flux max
p_p=p_prime(i,j); 

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';


%west south
i=ny;
j=1;
 
%get central and neighborhood nodes;
p_w=0; %Inlet
p_n=p_prime(i-1,j);
p_e=p_prime(i,j+1);
p_s=0;% No perpendicular flux max
p_p=p_prime(i,j); 

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p_prime(i,j,:,:)=(wlsqOp*difp_vec)';
