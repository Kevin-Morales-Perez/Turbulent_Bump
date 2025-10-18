function [grad_tau_xx] = computeTauGradient(tau_xx,wlsqOperator,...
    grad_tau_xx)
%Computation of gradients at the cell centers

%Turbulent Stresses
size_fields= size(tau_xx);
ny=size_fields(1);
nx=size_fields(2);


%internal cells 
for i= 2:ny-1
    for j=2:nx-1

        %get central and neighborhood nodes;
        tau_xx_w=tau_xx(i,j-1);
        tau_xx_n=tau_xx(i-1,j);
        tau_xx_e=tau_xx(i,j+1);
        tau_xx_s=tau_xx(i+1,j);
        tau_xx_p=tau_xx(i,j); 

        diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau
        % difference vec

        %Get gradient operator
        wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

        %Compute gradient
        grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';
       
    end
end

% IN EDGES AND CORNERS CONSIDER BOUNDARY CONDITIONS FOR tau_mn!  _____

%________________________EDGES_________________________
%West Edge (
j=1;
for i=2:ny-1

    %get central and neighborhood nodes;

    tau_xx_w=0; %Inlet
    tau_xx_n=tau_xx(i-1,j);
    tau_xx_e=tau_xx(i,j+1);
    tau_xx_s=tau_xx(i+1,j);
    tau_xx_p=tau_xx(i,j);

    diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau      
    % difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';

end


%North Edge
i=1;
for j=2:nx-1

    %get central and neighborhood nodes;
    tau_xx_w=tau_xx(i,j-1);
    tau_xx_n=tau_xx(i,j); %Neumman
    tau_xx_e=tau_xx(i,j+1);
    tau_xx_s=tau_xx(i+1,j);
    tau_xx_p=tau_xx(i,j); 

    diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau     
    %     % difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';
    
end


%East Edge
j=nx;
for i=2:ny-1
    
    %get central and neighborhood nodes;
    tau_xx_w=tau_xx(i,j-1);
    tau_xx_n=tau_xx(i-1,j);
    tau_xx_e=tau_xx(i,j);%outlet condition (Neumman)
    tau_xx_s=tau_xx(i+1,j);
    tau_xx_p=tau_xx(i,j); 

    diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau  
    %        % difference vec

    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

    %Compute gradient
    grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';


end


%South Edge
i=ny;
for j=2:nx-1
    
    %over the plate
    %get central and neighborhood nodes;
    tau_xx_w=tau_xx(i,j-1);
    tau_xx_n=tau_xx(i-1,j);
    tau_xx_e=tau_xx(i,j+1);
    tau_xx_s=0;% No turbulent stress 
    tau_xx_p=tau_xx(i,j); 
    
    diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau        
    %  % difference vec
    
    %Get gradient operator
    wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);
    
    %Compute gradient
    grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';

end


%______________________Corners______________________________
%west north
i=1;
j=1;

%get central and neighborhood nodes;
tau_xx_w=0; %Inlet
tau_xx_n=tau_xx(i,j); %Neumman
tau_xx_e=tau_xx(i,j+1);
tau_xx_s=tau_xx(i+1,j);
tau_xx_p=tau_xx(i,j); 

diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau 
%         % difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';


%east north
i=1;
j=nx;

%get central and neighborhood nodes;

tau_xx_w=tau_xx(i,j-1);
tau_xx_n=tau_xx(i,j); %Neumman
tau_xx_e=tau_xx(i,j);%outlet condition
tau_xx_s=tau_xx(i+1,j);
tau_xx_p=tau_xx(i,j);

diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau         % difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';


%east soth
i=ny;
j=nx;

%get central and neighborhood nodes;
tau_xx_w=tau_xx(i,j-1);
tau_xx_n=tau_xx(i-1,j);
tau_xx_e=tau_xx(i,j);%outlet condition
tau_xx_s=0;% Symetric condition
tau_xx_p=tau_xx(i,j); 

diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau         % difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';


%west south
i=ny;
j=1;
 
%get central and neighborhood nodes;
tau_xx_w=0; %Inlet
tau_xx_n=tau_xx(i-1,j);
tau_xx_e=tau_xx(i,j+1);
tau_xx_s=tau_xx(i,j);% Symetric condition
tau_xx_p=tau_xx(i,j); 

diftau_vec=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;%tau         % difference vec

%Get gradient operator
wlsqOp =reshape(wlsqOperator(i,j,:,:),[2,4]);

%Compute gradient
grad_tau_xx(i,j,:,:)=(wlsqOp*diftau_vec)';