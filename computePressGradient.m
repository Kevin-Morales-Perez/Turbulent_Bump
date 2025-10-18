function [grad_p] = computePressGradient(p,wlsqOperator,grad_p,p0)
%Computation of gradients at the cell centers

%pressure 
size_fields= size(p);
ny=size_fields(1);
nx=size_fields(2);


%internal cells 
for i= 2:ny-1
    for j=2:nx-1

        %get central and neighborhood nodes;
        p_w=p(i,j-1);
        p_n=p(i-1,j);
        p_e=p(i,j+1);
        p_s=p(i+1,j);
        p_p=p(i,j); 

        difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

        %Compute gradient
        grad_p(i,j,:,:)=(wlsqOp*difp_vec)';
       
    end
end

% IN EDGES AND CORNERS CONSIDER BOUNDARY CONDITIONS FOR PRESSURE!  _____

%________________________EDGES_________________________
%West Edge (
j=1;
for i=2:ny-1

    %get central and neighborhood nodes;

    p_w=(3/2)*p(i,j) - 0.5*p(i,j+1); %Inlet
    p_n=p(i-1,j);
    p_e=p(i,j+1);
    p_s=p(i+1,j);
    p_p=p(i,j);

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p(i,j,:,:)=(wlsqOp*difp_vec)';


end


%North Edge
i=1;
for j=2:nx-1

    %get central and neighborhood nodes;
    p_w=p(i,j-1);
    p_n=p0;%(3/2)*p(i,j) - 0.5*p(i+1,j);%Neumman condition
    p_e=p(i,j+1);
    p_s=p(i+1,j);
    p_p=p(i,j); 

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p(i,j,:,:)=(wlsqOp*difp_vec)';
    
end


%East Edge
j=nx;
for i=2:ny-1
    
    %get central and neighborhood nodes;
    p_w=p(i,j-1);
    p_n=p(i-1,j);
    p_e=p0;%outlet condition
    p_s=p(i+1,j);
    p_p=p(i,j); 

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p(i,j,:,:)=(wlsqOp*difp_vec)';


end


%South Edge
i=ny;
for j=2:nx-1

    %get central and neighborhood nodes;
    p_w=p(i,j-1);
    p_n=p(i-1,j);
    p_e=p(i,j+1);
    p_s=p(i,j);% No perpendicular flux max
    p_p=p(i,j); 

    difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_p(i,j,:,:)=(wlsqOp*difp_vec)';

end


%______________________Corners______________________________
%west north
i=1;
j=1;

%get central and neighborhood nodes;
p_w=(3/2)*p(i,j) - 0.5*p(i,j+1); %Inlet
p_n=p0; %(3/2)*p(i,j) - 0.5*p(i+1,j);%Dirichlet condition (Neumman)
p_e=p(i,j+1);
p_s=p(i+1,j);
p_p=p(i,j); 

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p(i,j,:,:)=(wlsqOp*difp_vec)';


%east north
i=1;
j=nx;

%get central and neighborhood nodes;

p_w=p(i,j-1);
p_n=p0;%(3/2)*p(i,j) - 0.5*p(i+1,j);%Dirichlet condition (Neumman)
p_e=p0;%outlet condition
p_s=p(i+1,j);
p_p=p(i,j);

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p(i,j,:,:)=(wlsqOp*difp_vec)';


%east soth
i=ny;
j=nx;

%get central and neighborhood nodes;
p_w=p(i,j-1);
p_n=p(i-1,j);
p_e=p0;%outlet condition
p_s=p(i,j);% No perpendicular flux max
p_p=p(i,j); 

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p(i,j,:,:)=(wlsqOp*difp_vec)';


%west south
i=ny;
j=1;
 
%get central and neighborhood nodes;
p_w=(3/2)*p(i,j) - 0.5*p(i,j+1); %Inlet
p_n=p(i-1,j);
p_e=p(i,j+1);
p_s=p(i,j);% No perpendicular flux max
p_p=p(i,j); 

difp_vec=[p_w;p_n;p_e;p_s]-p_p;%pressure difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_p(i,j,:,:)=(wlsqOp*difp_vec)';
