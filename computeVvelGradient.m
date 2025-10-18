function [grad_v] = computeVvelGradient(v,wlsqOperator,grad_v)
%Computation of gradients at the cell centers

%pressure 
size_fields= size(v);
ny=size_fields(1);
nx=size_fields(2);


%internal cells 
for i= 2:ny-1
    for j=2:nx-1

        %get central and neighborhood nodes;
        v_w=v(i,j-1);
        v_n=v(i-1,j);
        v_e=v(i,j+1);
        v_s=v(i+1,j);
        v_p=v(i,j); 

        difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

        %Compute gradient
        grad_v(i,j,:,:)=(wlsqOp*difv_vec)';
       
    end
end

% IN EDGES AND CORNERS CONSIDER BOUNDARY CONDITIONS FOR U velocity!  _____

%________________________EDGES_________________________
%West Edge (
j=1;
for i=2:ny-1

    %get central and neighborhood nodes;

    v_w=0; %Inlet
    v_n=v(i-1,j);
    v_e=v(i,j+1);
    v_s=v(i+1,j);
    v_p=v(i,j);

    difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_v(i,j,:,:)=(wlsqOp*difv_vec)';


end


%North Edge
i=1;
for j=2:nx-1

    %get central and neighborhood nodes;
    v_w=v(i,j-1);
    v_n=v(i,j); %Neumman
    v_e=v(i,j+1);
    v_s=v(i+1,j);
    v_p=v(i,j); 

    difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_v(i,j,:,:)=(wlsqOp*difv_vec)';
    
end


%East Edge
j=nx;
for i=2:ny-1
    
    %get central and neighborhood nodes;
    v_w=v(i,j-1);
    v_n=v(i-1,j);
    v_e=v(i,j);%outlet condition (Neumman)
    v_s=v(i+1,j);
    v_p=v(i,j); 

    difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_v(i,j,:,:)=(wlsqOp*difv_vec)';


end


%South Edge
i=ny;
for j=2:nx-1
    
    %over the plate
    %get central and neighborhood nodes;
    v_w=v(i,j-1);
    v_n=v(i-1,j);
    v_e=v(i,j+1);
    v_s=0;% No mass flux 
    v_p=v(i,j); 
    
    difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec
    
    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);
    
    %Compute gradient
    grad_v(i,j,:,:)=(wlsqOp*difv_vec)';

end


%______________________Corners______________________________
%west north
i=1;
j=1;

%get central and neighborhood nodes;
v_w=0; %Inlet
v_n=v(i,j); %Neumman
v_e=v(i,j+1);
v_s=v(i+1,j);
v_p=v(i,j); 

difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_v(i,j,:,:)=(wlsqOp*difv_vec)';


%east north
i=1;
j=nx;

%get central and neighborhood nodes;

v_w=v(i,j-1);
v_n=v(i,j); %Neumman
v_e=v(i,j);%outlet condition
v_s=v(i+1,j);
v_p=v(i,j);

difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_v(i,j,:,:)=(wlsqOp*difv_vec)';


%east soth
i=ny;
j=nx;

%get central and neighborhood nodes;
v_w=v(i,j-1);
v_n=v(i-1,j);
v_e=v(i,j);%outlet condition
v_s=0;% Symetric condition
v_p=v(i,j); 

difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_v(i,j,:,:)=(wlsqOp*difv_vec)';


%west south
i=ny;
j=1;
 
%get central and neighborhood nodes;
v_w=0; %Inlet
v_n=v(i-1,j);
v_e=v(i,j+1);
v_s=v(i,j);% Symetric condition
v_p=v(i,j); 

difv_vec=[v_w;v_n;v_e;v_s]-v_p;%u vel difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_v(i,j,:,:)=(wlsqOp*difv_vec)';