function [p_prime,rsid_p,err_p] = solve_presscorr(err_p,max_iterations_p...
    ,nx,ny,p_prime,ap_P,ap_W,ap_N,ap_E,ap_S,suP,rsid_p...
    ,epsilon_p)
%solution of pressure correction equation
%   ApP*PP'=ApW*PW' + ApN*PN' + ApE*PE' + ApS*PS' + SuP

    %use for this Gauss - Seidel 
    
    iterations_p=0;
    
    while err_p > epsilon_p && iterations_p < max_iterations_p
        iterations_p=iterations_p+1;
    
        %Top wall (North)
        i=1;
        for j=2:nx-1
    
            p_prime(i,j)= (ap_W(i,j)*p_prime(i,j-1) + ...
                ap_E(i,j)*p_prime(i,j+1) + ap_S(i,j)*p_prime(i+1,j) + ...
                suP(i,j))*(1/ap_P(i,j));
        end
    
        %Right wall (East)
        j=nx;
        for i=2:ny-1
    
            p_prime(i,j)= (ap_W(i,j)*p_prime(i,j-1) + ...
                ap_N(i,j)*p_prime(i-1,j) + ap_S(i,j)*p_prime(i+1,j) + ...
                suP(i,j))*(1/ap_P(i,j));
        end
    
        %Bottom wall (South)
        i=ny;
        for j=2:nx-1
            p_prime(i,j)= (ap_W(i,j)*p_prime(i,j-1) + ... 
                ap_N(i,j)*p_prime(i-1,j) + ap_E(i,j)*p_prime(i,j+1) + ...
                suP(i,j))*(1/ap_P(i,j));
        end
    
        %Left wall (West)
    
        j=1;
        for i =2:ny-1
    
            p_prime(i,j)= (ap_N(i,j)*p_prime(i-1,j) + ...
                ap_E(i,j)*p_prime(i,j+1) + ap_S(i,j)*p_prime(i+1,j) + ... 
                suP(i,j))*(1/ap_P(i,j));
        end
    
        %Corners 
       
        
        %North - West
        i=1;
        j=1;
    
        p_prime(i,j)=(ap_E(i,j)*p_prime(i,j+1) + ... 
            ap_S(i,j)*p_prime(i+1,j) + suP(i,j))*(1/ap_P(i,j));
    
        %North - East
        i=1;
        j=nx;
    
        p_prime(i,j)= (ap_W(i,j)*p_prime(i,j-1) + ... 
            ap_S(i,j)*p_prime(i+1,j) + suP(i,j))*(1/ap_P(i,j));
        
    
        %South - East
        i=ny;
        j=nx;
    
        p_prime(i,j)= (ap_W(i,j)*p_prime(i,j-1) + ... 
            ap_N(i,j)*p_prime(i-1,j) + suP(i,j))*(1/ap_P(i,j));
        
        %South - West
        i=ny;
        j=1;
    
        p_prime(i,j)= (ap_N(i,j)*p_prime(i-1,j) + ... 
            ap_E(i,j)*p_prime(i,j+1) + suP(i,j))*(1/ap_P(i,j));
    
        %Interior cells
        for i=2:ny-1
            for j=2:nx-1
                p_prime(i,j)= (ap_W(i,j)*p_prime(i,j-1) + ...
                    ap_N(i,j)*p_prime(i-1,j) + ap_E(i,j)*p_prime(i,j+1) + ...
                    ap_S(i,j)*p_prime(i+1,j) + suP(i,j))*(1/ap_P(i,j));
            end
        end
    
        %Residual computation
        rsid_p=zeros(ny,nx);


       
    
        %Walls
    
        %West wall
        j=1;
        for i =2:ny-1
            rsid_p(i,j) = p_prime(i,j) - ...
                (ap_N(i,j)*p_prime(i-1,j) + ... 
                ap_E(i,j)*p_prime(i,j+1) + ap_S(i,j)*p_prime(i+1,j) + ...
                suP(i,j))*(1/ap_P(i,j));
        end
    
        %North wall
    
        i=1;
        for j=2:nx-1
            rsid_p(i,j) = p_prime(i,j) - ...
              (ap_W(i,j)*p_prime(i,j-1) + ... 
                ap_E(i,j)*p_prime(i,j+1) + ap_S(i,j)*p_prime(i+1,j) + ... 
                suP(i,j))*(1/ap_P(i,j));
        end
    
        %East wall
        j=nx;
        for i=2:ny-1
    
            rsid_p(i,j) = p_prime(i,j) - ... 
                (ap_W(i,j)*p_prime(i,j-1) + ... 
                ap_N(i,j)*p_prime(i-1,j) + ap_S(i,j)*p_prime(i+1,j) + ... 
                suP(i,j))*(1/ap_P(i,j));
    
        end
    
    
        %South wall
        i=ny;
        for j=2:nx-1
    
            rsid_p(i,j) = p_prime(i,j) - ...
               (ap_W(i,j)*p_prime(i,j-1) + ... 
                ap_N(i,j)*p_prime(i-1,j) + ap_E(i,j)*p_prime(i,j+1) + ... 
                suP(i,j))*(1/ap_P(i,j));
    
        end
    
        %Corners
        
        %North - West
        i=1;
        j=1;
    
        rsid_p(i,j) = p_prime(i,j) - ... 
            (ap_E(i,j)*p_prime(i,j+1) + ... 
            ap_S(i,j)*p_prime(i+1,j) + suP(i,j))*(1/ap_P(i,j));
    
        %North - East
        i=1;
        j=nx;
    
        rsid_p(i,j) = p_prime(i,j) - ... 
            (ap_W(i,j)*p_prime(i,j-1) + ... 
            ap_S(i,j)*p_prime(i+1,j) + suP(i,j))*(1/ap_P(i,j));
    
        %South - East
        i=ny;
        j=nx;
    
        rsid_p(i,j) = p_prime(i,j) - ... 
            (ap_W(i,j)*p_prime(i,j-1) + ... 
            ap_N(i,j)*p_prime(i-1,j) + suP(i,j))*(1/ap_P(i,j));
    
        %South - West
        i=ny;
        j=1;
    
        rsid_p(i,j) = p_prime(i,j) - ... 
            (ap_N(i,j)*p_prime(i-1,j) +  ... 
            ap_E(i,j)*p_prime(i,j+1) + suP(i,j))*(1/ap_P(i,j));
    
        %Interior cells 
    
        for i=2:ny-1
            for j=2:nx-1
    
                rsid_p(i,j) = p_prime(i,j) - ...
                    (ap_W(i,j)*p_prime(i,j-1) + ... 
                    ap_N(i,j)*p_prime(i-1,j) + ap_E(i,j)*p_prime(i,j+1) + ... 
                    ap_S(i,j)*p_prime(i+1,j) + suP(i,j))*(1/ap_P(i,j));
    
            end
        end
     
    
        %Error from residual
        err_p=rms(rsid_p(:));
        
    end
        
end