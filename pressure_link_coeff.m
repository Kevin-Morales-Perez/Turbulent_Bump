function [ap_W,ap_N,ap_E,ap_S,ap_P,suP] = pressure_link_coeff(u_face,...
    v_face,aP,aPv,ap_W,ap_N,ap_E,ap_S,ap_P,suP,weightDistFactors,...
    cellVols,distNeighbNods,lgtFaces,nx_upstr,nx_dwnstr)
%Pressure correction coeffitients for non -orthogonal mesh
%pressure correction coeffitients ap_W,ap_N,ap_E,ap_S,ap_P
%the main equation is AP_p*p_p'= AP_w*p_w' + AP_n*p_n' + AP_e*p_e' 
%+Ap_S*p_s' + suP
%each coeffitient is [wf*CelVolp/aP_p + (1- wf)*CelVol_a/aP_a]*(fAa/dist_a)
%suP: Mass inbalance
%u_face,v_face: Velocities normal to faces
%aP,aPv:Coeffitient from momentum equations
%weightDistFactors: Weight distance factors
%cellVols: volumes of the cells
%distNeighbNods: Distances between cell centers and neigborhood nodes
%lgtFaces: Area of the faces of each cell  (since it is a 2D code are lenghts)
 %nx_dwnstr: number of cells in the free downstream zone 

    %sizes
    size_field=size(aP);
    nx=size_field(2);
    ny=size_field(1);

    %Interior cells except cells near SBC

    for i =2:ny-2
        for j = 2:nx-1

            %get weight distance factors

            weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

            %split factors
            weight_w=weight_F(1);
            weight_n=weight_F(2);
            weight_e=weight_F(3);
            weight_s=weight_F(4);

            %get volumes of central and nb cells 
            volW=cellVols(i,j-1); 
            volN=cellVols(i-1,j);
            volE=cellVols(i,j+1);
            volS=cellVols(i+1,j);
            volP=cellVols(i,j);

            %get distances from central node to nb nodes
            
            distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

            %split distNods

            dist_w=distNods(1);
            dist_n=distNods(2);
            dist_e=distNods(3);
            dist_s=distNods(4);

            %get cell face areas

            faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

            %split areas

            fAw=faceAreas(1);
            fAn=faceAreas(2);
            fAe=faceAreas(3);
            fAs=faceAreas(4);

            % Calculate pressure correction coefficients for nb cells
            ap_W(i,j) =(fAw/dist_w)*((volP/aP(i,j))*weight_w  +...
                (1-weight_w)*(volW/aP(i,j-1)));

            ap_N(i,j) = (fAn/dist_n)*((volP/aP(i,j))*weight_n  + ...
                (1-weight_n)*(volN/aP(i-1,j)));

            ap_E(i,j) = (fAe/dist_e)*((volP/aP(i,j))*weight_e  + ...
                (1-weight_e)*(volE/aP(i,j+1)));

            ap_S(i,j) = (fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
                (1-weight_s)*(volS/aP(i+1,j)));

            %calculate pressure correction coeffitient for central cell

            ap_P(i,j)= ap_W(i,j) + ap_N(i,j) + ap_E(i,j) + ap_S(i,j);

            %calculate source for pressure correction (mass inbalance)

            suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
                fAn*v_face(i,j) -fAs*v_face(i+1,j));
        end
    end

    %Cell near SBC
    
    i=ny-1;
    for j=2:nx-1

        %get weight distance factors

        weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

        %split factors
        weight_w=weight_F(1);
        weight_n=weight_F(2);
        weight_e=weight_F(3);
        weight_s=weight_F(4);

        %get volumes of central and nb cells 
        volW=cellVols(i,j-1); 
        volN=cellVols(i-1,j);
        volE=cellVols(i,j+1);
        volS=cellVols(i+1,j);
        volP=cellVols(i,j);

        %get distances from central node to nb nodes
        
        distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

        %split distNods

        dist_w=distNods(1);
        dist_n=distNods(2);
        dist_e=distNods(3);
        dist_s=distNods(4);

        %get cell face areas

        faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

        %split areas

        fAw=faceAreas(1);
        fAn=faceAreas(2);
        fAe=faceAreas(3);
        fAs=faceAreas(4);


        if (j<=nx_upstr) || j>(nx-nx_dwnstr)
            %cells over the plate use aPv
            ap_S(i,j) = (fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
                (1-weight_s)*(volS/aPv(j)));

        else
            %free stream
            ap_S(i,j) =(fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
                (1-weight_s)*(volS/aP(i+1,j)));
        end

        % Calculate pressure correction coefficients for nb cells
        ap_W(i,j) = (fAw/dist_w)*((volP/aP(i,j))*weight_w  +...
            (1-weight_w)*(volW/aP(i,j-1)));

        ap_N(i,j) =  (fAn/dist_n)*((volP/aP(i,j))*weight_n  + ...
            (1-weight_n)*(volN/aP(i-1,j)));

        ap_E(i,j) = (fAe/dist_e)*((volP/aP(i,j))*weight_e  + ...
            (1-weight_e)*(volE/aP(i,j+1)));

        %calculate pressure correction coeffitient for central cell

        ap_P(i,j)= ap_W(i,j) + ap_N(i,j) + ap_E(i,j) + ap_S(i,j);

        %calculate source for pressure correction (mass inbalance)

        suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
            fAn*v_face(i,j) -fAs*v_face(i+1,j));

    end

   %________________________EDGES_________________________
    %West Edge
    j=1;
    for i=2:ny-1

        %get weight distance factors

        weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

        %split factors
  
        weight_n=weight_F(2);
        weight_e=weight_F(3);
        weight_s=weight_F(4);

        %get volumes of central and nb cells 
 
        volN=cellVols(i-1,j);
        volE=cellVols(i,j+1);
        volS=cellVols(i+1,j);
        volP=cellVols(i,j);

        %get distances from central node to nb nodes
        
        distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

        %split distNods

        dist_n=distNods(2);
        dist_e=distNods(3);
        dist_s=distNods(4);

        %get cell face areas

        faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

        %split areas

        fAw=faceAreas(1);
        fAn=faceAreas(2);
        fAe=faceAreas(3);
        fAs=faceAreas(4);

        % Calculate pressure correction coefficients for nb cells
       

        ap_N(i,j) = (fAn/dist_n)*((volP/aP(i,j))*weight_n  + ...
            (1-weight_n)*(volN/aP(i-1,j)));

        ap_E(i,j) = (fAe/dist_e)*((volP/aP(i,j))*weight_e  + ...
            (1-weight_e)*(volE/aP(i,j+1)));

        ap_S(i,j) = (fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
            (1-weight_s)*(volS/aP(i+1,j)));

        %calculate pressure correction coeffitient for central cell

        ap_P(i,j)= ap_N(i,j) + ap_E(i,j) + ap_S(i,j);

        %calculate source for pressure correction (mass inbalance)

        suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
            fAn*v_face(i,j) -fAs*v_face(i+1,j));

    end
    %North Edge
    i=1;
    for j=2:nx-1

        %get weight distance factors

        weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

        %split factors
        weight_w=weight_F(1);
        weight_e=weight_F(3);
        weight_s=weight_F(4);

        %get volumes of central and nb cells 
        volW=cellVols(i,j-1); 
        volE=cellVols(i,j+1);
        volS=cellVols(i+1,j);
        volP=cellVols(i,j);

        %get distances from central node to nb nodes
        
        distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

        %split distNods

        dist_w=distNods(1);
        dist_n=distNods(2);
        dist_e=distNods(3);
        dist_s=distNods(4);

        %get cell face areas

        faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

        %split areas

        fAw=faceAreas(1);
        fAn=faceAreas(2);
        fAe=faceAreas(3);
        fAs=faceAreas(4);

        % Calculate pressure correction coefficients for nb cells
        ap_W(i,j) = (fAw/dist_w)*((volP/aP(i,j))*weight_w  +...
            (1-weight_w)*(volW/aP(i,j-1)));

        ap_E(i,j) = (fAe/dist_e)*((volP/aP(i,j))*weight_e  + ...
            (1-weight_e)*(volE/aP(i,j+1)));

        ap_S(i,j) = (fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
            (1-weight_s)*(volS/aP(i+1,j))) -...
            (fAn/(4*dist_n))*(volP/aP(i,j));

        %calculate pressure correction coeffitient for central cell

        ap_P(i,j)= ap_W(i,j) + ap_E(i,j) + ap_S(i,j) + ...
            (fAn/(2*dist_n))*(volP/aP(i,j));

        %calculate source for pressure correction (mass inbalance)

        suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
            fAn*v_face(i,j) -fAs*v_face(i+1,j));

    end
    %East Edge
    j=nx;
    for i=2:ny-1

        %get weight distance factors

        weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

        %split factors
        weight_w=weight_F(1);
        weight_n=weight_F(2);
        %weight_e=weight_F(3);
        weight_s=weight_F(4);

        %get volumes of central and nb cells 
        volW=cellVols(i,j-1); 
        volN=cellVols(i-1,j);
        %volE=cellVols(i,j+1);
        volS=cellVols(i+1,j);
        volP=cellVols(i,j);

        %get distances from central node to nb nodes
        
        distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

        %split distNods

        dist_w=distNods(1);
        dist_n=distNods(2);
        dist_e=distNods(3);
        dist_s=distNods(4);

        %get cell face areas

        faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

        %split areas

        fAw=faceAreas(1);
        fAn=faceAreas(2);
        fAe=faceAreas(3);
        fAs=faceAreas(4);

        % Calculate pressure correction coefficients for nb cells
        ap_W(i,j) = (fAw/dist_w)*((volP/aP(i,j))*weight_w  +...
            (1-weight_w)*(volW/aP(i,j-1))) - ...
            (fAe/(4*dist_e))*(volP/aP(i,j));
        %corrected term :(fAe/(4*dist_e))*(volP/aP(i,j))
        
        %old term :volP/(4*aP(i,j)*dist_e)

        ap_N(i,j) = (fAn/dist_n)*((volP/aP(i,j))*weight_n  + ...
            (1-weight_n)*(volN/aP(i-1,j)));

        ap_S(i,j) =(fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
                (1-weight_s)*(volS/aP(i+1,j)));

        %calculate pressure correction coeffitient for central cell

        ap_P(i,j)= ap_W(i,j) + ap_N(i,j) + ap_E(i,j) + ap_S(i,j) ... 
            + (fAe/(2*dist_e))*(volP/aP(i,j));

        %calculate source for pressure correction (mass inbalance)

        suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
            fAn*v_face(i,j) -fAs*v_face(i+1,j));
    end


    %South Edge (use here aPV from the symetric boundary condition)

    i=ny;
    for j=2:nx-1

        %get weight distance factors

        weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

        %split factors
        weight_w=weight_F(1);
        weight_n=weight_F(2);
        weight_e=weight_F(3);

        %get volumes of central and nb cells 
        volW=cellVols(i,j-1); 
        volN=cellVols(i-1,j);
        volE=cellVols(i,j+1);
        volP=cellVols(i,j);

        %get distances from central node to nb nodes
        
        distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

        %split distNods

        dist_w=distNods(1);
        dist_n=distNods(2);
        dist_e=distNods(3);
 

        %get cell face areas

        faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

        %split areas

        fAw=faceAreas(1);
        fAn=faceAreas(2);
        fAe=faceAreas(3);
        fAs=faceAreas(4);

        % Calculate pressure correction coefficients for nb cells
        ap_W(i,j) = (fAw/dist_w)*((volP/aP(i,j))*weight_w  +...
            (1-weight_w)*(volW/aP(i,j-1)));

        if (j<=nx_upstr) || j>(nx-nx_dwnstr) %over the plate

            ap_N(i,j) = (fAn/dist_n)*((volP/aPv(j))*weight_n  + ...
                (1-weight_n)*(volN/aP(i-1,j)));
        else
            ap_N(i,j) = (fAn/dist_n)*((volP/aP(i,j))*weight_n  + ...
                (1-weight_n)*(volN/aP(i-1,j)));
        end
        ap_E(i,j) = (fAe/dist_e)*((volP/aP(i,j))*weight_e  + ...
                (1-weight_e)*(volE/aP(i,j+1)));

        %calculate pressure correction coeffitient for central cell

        ap_P(i,j)= ap_W(i,j) + ap_N(i,j) + ap_E(i,j);

        %calculate source for pressure correction (mass inbalance)

        suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
            fAn*v_face(i,j) -fAs*v_face(i+1,j));

    end
    %______________________Corners______________________________
    %west north
    i=1;
    j=1;

    %get weight distance factors

    weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

    %split factors
    weight_e=weight_F(3);
    weight_s=weight_F(4);

    %get volumes of central and nb cells 
    volE=cellVols(i,j+1);
    volS=cellVols(i+1,j);
    volP=cellVols(i,j);

    %get distances from central node to nb nodes
    
    distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

    %split distNods

    dist_n=distNods(2);
    dist_e=distNods(3);
    dist_s=distNods(4);

    %get cell face areas

    faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

    %split areas

    fAw=faceAreas(1);
    fAn=faceAreas(2);
    fAe=faceAreas(3);
    fAs=faceAreas(4);

    % Calculate pressure correction coefficients for nb cells

    ap_E(i,j) =(fAe/dist_e)*((volP/aP(i,j))*weight_e  + ...
            (1-weight_e)*(volE/aP(i,j+1)));

    ap_S(i,j) = (fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
            (1-weight_s)*(volS/aP(i+1,j))) - ...
            (fAn/(4*dist_n))*(volP/aP(i,j));

    %calculate pressure correction coeffitient for central cell

    ap_P(i,j)= ap_E(i,j) + ap_S(i,j) + (fAn/(2*dist_n))*(volP/aP(i,j));

    %calculate source for pressure correction (mass inbalance)

    suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
        fAn*v_face(i,j) -fAs*v_face(i+1,j));


    %east north
    i=1;
    j=nx;

    %get weight distance factors

    weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

    %split factors
    weight_w=weight_F(1);
    weight_s=weight_F(4);

    %get volumes of central and nb cells 
    volW=cellVols(i,j-1); 
    volS=cellVols(i+1,j);
    volP=cellVols(i,j);

    %get distances from central node to nb nodes
    
    distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

    %split distNods

    dist_w=distNods(1);
    dist_e=distNods(3);
    dist_s=distNods(4);

    %get cell face areas

    faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

    %split areas

    fAw=faceAreas(1);
    fAn=faceAreas(2);
    fAe=faceAreas(3);
    fAs=faceAreas(4);

    % Calculate pressure correction coefficients for nb cells
    ap_W(i,j) = (fAw/dist_w)*((volP/aP(i,j))*weight_w  +...
        (1-weight_w)*(volW/aP(i,j-1))) -...
        (fAe/(4*dist_e))*(volP/aP(i,j));

    ap_S(i,j) = (fAs/dist_s)*((volP/aP(i,j))*weight_s  + ... 
            (1-weight_s)*(volS/aP(i+1,j))) - ...
            (fAn/(4*dist_n))*(volP/aP(i,j));

    %calculate pressure correction coeffitient for central cell

    ap_P(i,j)= ap_W(i,j) + ap_S(i,j) + (fAe/(2*dist_e))*(volP/aP(i,j)) ...
        + (fAn/(2*dist_n))*(volP/aP(i,j));

    %calculate source for pressure correction (mass inbalance)

    suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
        fAn*v_face(i,j) -fAs*v_face(i+1,j));

    %east south
    i=ny;
    j=nx;

    %get weight distance factors

    weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

    %split factors
    weight_w=weight_F(1);
    weight_n=weight_F(2);

    %get volumes of central and nb cells 
    volW=cellVols(i,j-1); 
    volN=cellVols(i-1,j);
    volP=cellVols(i,j);

    %get distances from central node to nb nodes
    
    distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

    %split distNods

    dist_w=distNods(1);
    dist_n=distNods(2);
    dist_e=distNods(3);

    %get cell face areas

    faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

    %split areas

    fAw=faceAreas(1);
    fAn=faceAreas(2);
    fAe=faceAreas(3);
    fAs=faceAreas(4);

    % Calculate pressure correction coefficients for nb cells
    ap_W(i,j) =(fAw/dist_w)*((volP/aP(i,j))*weight_w  +...
            (1-weight_w)*(volW/aP(i,j-1))) - ...
            (fAe/(4*dist_e))*(volP/aP(i,j));

    ap_N(i,j) = (fAn/dist_n)*((volP/aPv(j))*weight_n  + ...
            (1-weight_n)*(volN/aP(i-1,j)));

    %calculate pressure correction coeffitient for central cell

    ap_P(i,j)= ap_W(i,j) + ap_N(i,j) + (fAe/(2*dist_e))*(volP/aP(i,j));

    %calculate source for pressure correction (mass inbalance)

    suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
        fAn*v_face(i,j) -fAs*v_face(i+1,j));


    %west south
    i=ny;
    j=1;

    %get weight distance factors

    weight_F=reshape(weightDistFactors(i,j,:),[1,4]);

    %split factors
    weight_n=weight_F(2);
    weight_e=weight_F(3);


    %get volumes of central and nb cells  
    volN=cellVols(i-1,j);
    volE=cellVols(i,j+1);
    volP=cellVols(i,j);

    %get distances from central node to nb nodes
    
    distNods=reshape(distNeighbNods(i,j,:,:),[1,4]);

    %split distNods

    dist_n=distNods(2);
    dist_e=distNods(3);

    %get cell face areas

    faceAreas=reshape(lgtFaces(i,j,:),[1 4]);

    %split areas

    fAw=faceAreas(1);
    fAn=faceAreas(2);
    fAe=faceAreas(3);
    fAs=faceAreas(4);

    % Calculate pressure correction coefficients for nb cells
    
    ap_N(i,j) = (fAn/dist_n)*((volP/aPv(j))*weight_n  + ...
            (1-weight_n)*(volN/aP(i-1,j)));

    ap_E(i,j) = (fAe/dist_e)*((volP/aP(i,j))*weight_e  + ...
            (1-weight_e)*(volE/aP(i,j+1)));

    %calculate pressure correction coeffitient for central cell

    ap_P(i,j)= ap_N(i,j) + ap_E(i,j);

    %calculate source for pressure correction (mass inbalance)

    suP(i,j) = -(-fAw*u_face(i,j) + fAe*u_face(i,j+1) + ...
        fAn*v_face(i,j) -fAs*v_face(i+1,j));

end