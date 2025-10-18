function [u,rsid_x,err_x] = solve_momentum(err_x,max_iterations_uv,...
    xGridPoints,yGridPoints,u,alpha_uv,aP,aW,aN,aE,aS,suX,u_star,rsid_x)
%Solution of Momentum equations 
%   Solution of matrix momentum equation using Gauss - Seidel Method with
%underelaxation to obtain velocities,the systems are in the form of: 
% Ap*phi_p = Aw*phi_w + An*phi_n + Ae*phi_e + As*phi_s + Su, where phi is 
%the transported variable
    
    %Residual
    
    iterations_x=0;
    
    while iterations_x < max_iterations_uv %err_x > epsilon_uv && iterations_x < max_iterations_uv
        iterations_x=iterations_x+1;
    
        %Top wall (North)
        i=1;
        for j=2:xGridPoints-1
    
            u(i,j)= (alpha_uv*(aW(i,j)*u(i,j-1) + aE(i,j)*u(i,j+1) + ... 
                aS(i,j)*u(i+1,j) + suX(i,j))*(1/aP(i,j))) + ...
                (1-alpha_uv)*u_star(i,j);
        end
    
        %Right wall (East)
        j=xGridPoints;
        for i=2:yGridPoints-1
    
            u(i,j)= (alpha_uv*(aW(i,j)*u(i,j-1) + aN(i,j)*u(i-1,j) + ... 
                aS(i,j)*u(i+1,j) + suX(i,j))*(1/aP(i,j))) ... 
                + (1-alpha_uv)*u_star(i,j);
        end
    
        %Bottom wall (South)
        i=yGridPoints;
        for j=2:xGridPoints-1
            u(i,j)= (alpha_uv*(aW(i,j)*u(i,j-1) + aN(i,j)*u(i-1,j) + ...
                aE(i,j)*u(i,j+1) + suX(i,j))*(1/aP(i,j))) + ... 
                (1-alpha_uv)*u_star(i,j);
        end
    
        %Left wall (West)
    
        j=1;
        for i =2:yGridPoints-1
    
            u(i,j)= (alpha_uv*(aN(i,j)*u(i-1,j) + aE(i,j)*u(i,j+1) + ... 
                aS(i,j)*u(i+1,j) + suX(i,j))*(1/aP(i,j))) + ... 
                (1-alpha_uv)*u_star(i,j);
        end
    
        %Corners 
        
        %North - West
        i=1;
        j=1;
    
        u(i,j)= (alpha_uv*(aE(i,j)*u(i,j+1) + aS(i,j)*u(i+1,j) + ... 
                    suX(i,j))*(1/aP(i,j))) + (1-alpha_uv)*u_star(i,j);
    
        %North - East
        i=1;
        j=xGridPoints;
    
        u(i,j)= (alpha_uv*(aW(i,j)*u(i,j-1) + aS(i,j)*u(i+1,j) + ... 
                suX(i,j))*(1/aP(i,j))) + (1-alpha_uv)*u_star(i,j);
    
        %South - East
        i=yGridPoints;
        j=xGridPoints;
    
        u(i,j)= (alpha_uv*(aW(i,j)*u(i,j-1) + aN(i,j)*u(i-1,j) + ...
                suX(i,j))*(1/aP(i,j))) + (1-alpha_uv)*u_star(i,j);
        
        %South - West
        i=yGridPoints;
        j=1;
    
        u(i,j)= (alpha_uv*(aN(i,j)*u(i-1,j) + ...
            aE(i,j)*u(i,j+1) + suX(i,j))*(1/aP(i,j))) + ... 
            (1-alpha_uv)*u_star(i,j);
    
        %Interior cells
        for i=2:yGridPoints-1
            for j=2:xGridPoints-1
                u(i,j)= (alpha_uv*(aW(i,j)*u(i,j-1) + aN(i,j)*u(i-1,j)+ ...
                    aE(i,j)*u(i,j+1) + aS(i,j)*u(i+1,j) + ... 
                    suX(i,j))*(1/aP(i,j))) + (1-alpha_uv)*u_star(i,j);
            end
        end
    
        %Raw Residual computation

        %Re-initialize residual after every iteration
        rsid_x=zeros(yGridPoints,xGridPoints);
    
        %Walls
    
        %West wall
        j=1;
        for i =2:yGridPoints-1
            rsid_x(i,j) = u(i,j) -(aN(i,j)*u(i-1,j) + ... 
                aE(i,j)*u(i,j+1) + aS(i,j)*u(i+1,j) + ...
                suX(i,j))*(1/aP(i,j));
        end
    
        %North wall
    
        i=1;
        for j=2:xGridPoints-1
            rsid_x(i,j) = u(i,j) -(aW(i,j)*u(i,j-1) + ... 
                aE(i,j)*u(i,j+1) + aS(i,j)*u(i+1,j) + ... 
                suX(i,j))*(1/aP(i,j));
        end
    
        %East wall
        j=xGridPoints;
        for i=2:yGridPoints-1
    
            rsid_x(i,j) = u(i,j) -(aW(i,j)*u(i,j-1) + ... 
                aN(i,j)*u(i-1,j) + aS(i,j)*u(i+1,j) + ... 
                suX(i,j))*(1/aP(i,j));
    
        end
    
    
        %South wall
        i=yGridPoints;
        for j=2:xGridPoints-1
    
            rsid_x(i,j) = u(i,j) -(aW(i,j)*u(i,j-1) + ... 
                aN(i,j)*u(i-1,j) + aE(i,j)*u(i,j+1) + ... 
                suX(i,j))*(1/aP(i,j));
    
        end
    
        %Corners
        
        %North - West
        i=1;
        j=1;
    
        rsid_x(i,j) = u(i,j) -(aE(i,j)*u(i,j+1) + ... 
            aS(i,j)*u(i+1,j) + suX(i,j))*(1/aP(i,j));
    
        %North - East
        i=1;
        j=xGridPoints;
    
        rsid_x(i,j) = u(i,j) -(aW(i,j)*u(i,j-1) + ... 
            aS(i,j)*u(i+1,j) + suX(i,j))*(1/aP(i,j));
    
        %South - East
        i=yGridPoints;
        j=xGridPoints;
    
        rsid_x(i,j) = u(i,j) -(aW(i,j)*u(i,j-1) + ... 
                aN(i,j)*u(i-1,j) + suX(i,j))*(1/aP(i,j));
    
        %South - West
        i=yGridPoints;
        j=1;
    
        rsid_x(i,j) = u(i,j) - (aN(i,j)*u(i-1,j) +  ... 
            aE(i,j)*u(i,j+1) + suX(i,j))*(1/aP(i,j));
    
        %Interior cells 
    
        for i=2:yGridPoints-1
            for j=2:xGridPoints-1
    
                rsid_x(i,j) = u(i,j) -(aW(i,j)*u(i,j-1) + ... 
                    aN(i,j)*u(i-1,j) + aE(i,j)*u(i,j+1) + aS(i,j)*u(i+1,j) +...
                    suX(i,j))*(1/aP(i,j));
    
            end
        end
    
        %Error from residual (RMS to avoid mesh size dependency)
        err_x=rms(rsid_x(:));
        
    end
    
end
