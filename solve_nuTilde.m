function [nu_tilde,rsid_nt,err_nt] = solve_nuTilde(err_nt,max_iterations_nt...
    ,nx,ny,nu_tilde,ant_P,ant_W,ant_N,ant_E,ant_S,suNT,rsid_nt...
    ,epsilon_nt)
%solution ofnu tilde linearized (Spallart Allmaras Equation)
%   ApP*PP'=ApW*PW' + ApN*PN' + ApE*PE' + ApS*PS' + SuP

    %use for this Gauss - Seidel 
    
    iterations_nt=0;
    
    while err_nt > epsilon_nt && iterations_nt < max_iterations_nt
        iterations_nt=iterations_nt+1;
    
        %Top wall (North)
        i=1;
        for j=2:nx-1
    
            nu_tilde(i,j)= (ant_W(i,j)*nu_tilde(i,j-1) + ...
                ant_E(i,j)*nu_tilde(i,j+1) + ant_S(i,j)*nu_tilde(i+1,j) + ...
                suNT(i,j))*(1/ant_P(i,j));
        end
    
        %Right wall (East)
        j=nx;
        for i=2:ny-1
    
            nu_tilde(i,j)= (ant_W(i,j)*nu_tilde(i,j-1) + ...
                ant_N(i,j)*nu_tilde(i-1,j) + ant_S(i,j)*nu_tilde(i+1,j) + ...
                suNT(i,j))*(1/ant_P(i,j));
        end
    
        %Bottom wall (South)
        i=ny;
        for j=2:nx-1
            nu_tilde(i,j)= (ant_W(i,j)*nu_tilde(i,j-1) + ... 
                ant_N(i,j)*nu_tilde(i-1,j) + ant_E(i,j)*nu_tilde(i,j+1) + ...
                suNT(i,j))*(1/ant_P(i,j));
        end
    
        %Left wall (West)
    
        j=1;
        for i =2:ny-1
    
            nu_tilde(i,j)= (ant_N(i,j)*nu_tilde(i-1,j) + ...
                ant_E(i,j)*nu_tilde(i,j+1) + ant_S(i,j)*nu_tilde(i+1,j) + ... 
                suNT(i,j))*(1/ant_P(i,j));
        end
    
        %Corners 
       
        
        %North - West
        i=1;
        j=1;
    
        nu_tilde(i,j)=(ant_E(i,j)*nu_tilde(i,j+1) + ... 
            ant_S(i,j)*nu_tilde(i+1,j) + suNT(i,j))*(1/ant_P(i,j));
    
        %North - East
        i=1;
        j=nx;
    
        nu_tilde(i,j)= (ant_W(i,j)*nu_tilde(i,j-1) + ... 
            ant_S(i,j)*nu_tilde(i+1,j) + suNT(i,j))*(1/ant_P(i,j));
        
    
        %South - East
        i=ny;
        j=nx;
    
        nu_tilde(i,j)= (ant_W(i,j)*nu_tilde(i,j-1) + ... 
            ant_N(i,j)*nu_tilde(i-1,j) + suNT(i,j))*(1/ant_P(i,j));
        
        %South - West
        i=ny;
        j=1;
    
        nu_tilde(i,j)= (ant_N(i,j)*nu_tilde(i-1,j) + ... 
            ant_E(i,j)*nu_tilde(i,j+1) + suNT(i,j))*(1/ant_P(i,j));
    
        %Interior cells
        for i=2:ny-1
            for j=2:nx-1
                nu_tilde(i,j)= (ant_W(i,j)*nu_tilde(i,j-1) + ...
                    ant_N(i,j)*nu_tilde(i-1,j) + ant_E(i,j)*nu_tilde(i,j+1) + ...
                    ant_S(i,j)*nu_tilde(i+1,j) + suNT(i,j))*(1/ant_P(i,j));
            end
        end
    
        %Residual computation
        rsid_nt=zeros(ny,nx);


       
    
        %Walls
    
        %West wall
        j=1;
        for i =2:ny-1
            rsid_nt(i,j) = nu_tilde(i,j) - ...
                (ant_N(i,j)*nu_tilde(i-1,j) + ... 
                ant_E(i,j)*nu_tilde(i,j+1) + ant_S(i,j)*nu_tilde(i+1,j) + ...
                suNT(i,j))*(1/ant_P(i,j));
        end
    
        %North wall
    
        i=1;
        for j=2:nx-1
            rsid_nt(i,j) = nu_tilde(i,j) - ...
              (ant_W(i,j)*nu_tilde(i,j-1) + ... 
                ant_E(i,j)*nu_tilde(i,j+1) + ant_S(i,j)*nu_tilde(i+1,j) + ... 
                suNT(i,j))*(1/ant_P(i,j));
        end
    
        %East wall
        j=nx;
        for i=2:ny-1
    
            rsid_nt(i,j) = nu_tilde(i,j) - ... 
                (ant_W(i,j)*nu_tilde(i,j-1) + ... 
                ant_N(i,j)*nu_tilde(i-1,j) + ant_S(i,j)*nu_tilde(i+1,j) + ... 
                suNT(i,j))*(1/ant_P(i,j));
    
        end
    
    
        %South wall
        i=ny;
        for j=2:nx-1
    
            rsid_nt(i,j) = nu_tilde(i,j) - ...
               (ant_W(i,j)*nu_tilde(i,j-1) + ... 
                ant_N(i,j)*nu_tilde(i-1,j) + ant_E(i,j)*nu_tilde(i,j+1) + ... 
                suNT(i,j))*(1/ant_P(i,j));
    
        end
    
        %Corners
        
        %North - West
        i=1;
        j=1;
    
        rsid_nt(i,j) = nu_tilde(i,j) - ... 
            (ant_E(i,j)*nu_tilde(i,j+1) + ... 
            ant_S(i,j)*nu_tilde(i+1,j) + suNT(i,j))*(1/ant_P(i,j));
    
        %North - East
        i=1;
        j=nx;
    
        rsid_nt(i,j) = nu_tilde(i,j) - ... 
            (ant_W(i,j)*nu_tilde(i,j-1) + ... 
            ant_S(i,j)*nu_tilde(i+1,j) + suNT(i,j))*(1/ant_P(i,j));
    
        %South - East
        i=ny;
        j=nx;
    
        rsid_nt(i,j) = nu_tilde(i,j) - ... 
            (ant_W(i,j)*nu_tilde(i,j-1) + ... 
            ant_N(i,j)*nu_tilde(i-1,j) + suNT(i,j))*(1/ant_P(i,j));
    
        %South - West
        i=ny;
        j=1;
    
        rsid_nt(i,j) = nu_tilde(i,j) - ... 
            (ant_N(i,j)*nu_tilde(i-1,j) +  ... 
            ant_E(i,j)*nu_tilde(i,j+1) + suNT(i,j))*(1/ant_P(i,j));
    
        %Interior cells 
    
        for i=2:ny-1
            for j=2:nx-1
    
                rsid_nt(i,j) = nu_tilde(i,j) - ...
                    (ant_W(i,j)*nu_tilde(i,j-1) + ... 
                    ant_N(i,j)*nu_tilde(i-1,j) + ant_E(i,j)*nu_tilde(i,j+1) + ... 
                    ant_S(i,j)*nu_tilde(i+1,j) + suNT(i,j))*(1/ant_P(i,j));
    
            end
        end
     
    
        %Error from residual
        err_nt=rms(rsid_nt(:));
        
    end