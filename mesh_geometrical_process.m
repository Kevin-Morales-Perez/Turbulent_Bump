function [cellVertxs,cellCentrs,lgtFaces,uVecParlFaces,faceCentrs,...
    cellVols,uVecNormFaces,distNeighbNods,uVecNeighbNods,...
    wlsqOperator,VecCentNodFaceCents,VecCentVertx,weightDistFactors,...
    Xctrs,Yctrs,distMinWall] =mesh_geometrical_process(X,Y,cellVertxs,cellCentrs,...
    lgtFaces,uVecParlFaces,faceCentrs,cellVols,uVecNormFaces,...
    distNeighbNods,uVecNeighbNods,wlsqOperator,VecCentNodFaceCents,...
    VecCentVertx,weightDistFactors,Xctrs,Yctrs,distMinWall,distPlat)

    %Get size of the mesh 
    
    [ny,nx]=size(X);
    nx=nx-1;
    ny=ny-1;

    %Compute cell vertex
    
    for i=1:ny %Iterate over rows
        for j=1:nx %Iterate over colums

            %Compute cell vertex

            Vwn=[X(i,j),Y(i,j)];%vertex west north
            Ven=[X(i,j+1),Y(i,j+1)];%vertex east north
            Ves=[X(i+1,j+1),Y(i+1,j+1)];%vertex east south
            Vws=[X(i+1,j),Y(i+1,j)];%vertex west south

            cellVertxs(i,j,:,:)=[Vwn;Ven;Ves;Vws];

            %_____________LEGACY CODE FOR CELL CENTERS_________________
            %(cellGeom.m)
            %1(1)=vws,1(2)=Vwn,1(3)=Ven,1(4)=Ves
            x1=[X(i+1,j) X(i,j) X(i,j+1) X(i+1,j+1)];%Define x vertex 
            % coordinates
     
            y1=[Y(i+1,j) Y(i,j) Y(i,j+1) Y(i+1,j+1)];%Define y vertex 
            % coordinates

            %funtion to calculate the centroid of each cell
            %the cell is divided in 2 triangles and one rectangle, then is
            %used the formula xc=sum(xi*ai)/a_total
            
            dx=x1(4)-x1(1);%step in X
            dy1= y1(2)-y1(1);%lenght of left face
            dy2=y1(3)-y1(4);%lenght of right face
            if y1(3)>y1(2)
                a=y1(4)-y1(1);%height of triangle 1
                b=y1(3)-y1(2);%height of triangle 2
                c= y1(2)-y1(4);% height of rectangle 1 
            
                xt1=(1/3)*dx; % x position of centroid for triangle 1
                xt2=(2/3)*dx; % x position of centroid for triangle 2
                y_in=y1(1);%y origin
    
            else
                a=y1(1)-y1(4);%height of triangle 1
                b=y1(2)-y1(3);%height of triangle 2
                c=y1(3)-y1(1);% height of rectangle 1 
            
                xt1=(2/3)*dx; % x position of centroid for triangle 1
                xt2=(1/3)*dx; % x position of centroid for triangle 2
                y_in=y1(4);%y origin
            end

            xt3=0.5*dx;% x position of centroid for rectangle 2 
            
            yt1=(2/3)*a;
            yt2=a + c + (1/3)*b;
            yt3= a + 0.5*c;
            
            at1=0.5*a*dx;%area of triangle 1;
            at2=0.5*b*dx;%area of triangle 2;
            at3=dx*c;%area of rectangle 3;
            
            a_total=0.5*(dy1+dy2)*dx;%total area of the cell
            xc=x1(1)+(xt1*at1 + xt2*at2 + xt3*at3)/a_total;%x centroid of the cell
            yc=y_in+ (at1*yt1 + at2*yt2 + at3*yt3)/a_total;%y centroid of the 

            cell_centroid=[xc,yc];
            
            cellCentrs(i,j,:)=cell_centroid;

            %___________________________________________________________

            %Use the cell vertexs to create vectors parallel to faces
            %VECTORS FOR FACE W AND E ARE ALLWAYS POSITIVE AND
            %ONLY WITH Y COMPONEN
            %VECTORS FOR  FACES N AND S HAVE ALWAYS POSITIVE X COMPONENT

            par_w=Vwn-Vws;
            par_n=Ven-Vwn;
            par_e=Ven-Ves;
            par_s=Ves-Vws;

            %calculate the length of each face

            lgt_w=vecnorm(par_w);
            lgt_n=vecnorm(par_n);
            lgt_e=vecnorm(par_e);
            lgt_s=vecnorm(par_s);
            %Store variables
            lgtFaces(i,j,:) = [lgt_w, lgt_n, lgt_e, lgt_s];

            %make unitary previous vectors;
            upar_w=par_w/lgt_w;
            upar_n=par_n/lgt_n;
            upar_e=par_w/lgt_e;
            upar_s=par_s/lgt_s;

            %Store unitary vectors paralel to faces
            uVecParlFaces(i,j,:,:)=[upar_w;upar_n;upar_e;upar_s];

            %Center of faces 

            %face w, starting point:vertex vws, 
            
            fcent_w=Vws +0.5*lgt_w*upar_w;
            fcent_n=Vwn +0.5*lgt_n*upar_n;
            fcent_e=Ves +0.5*lgt_e*upar_e;
            fcent_s=Vws +0.5*lgt_s*upar_s;

            faceCentrs(i,j,:,:)=[fcent_w;fcent_n;fcent_e;fcent_s];

            %Cell Vols
            %use trapezium formula

            cellVols(i,j)=(lgt_w+lgt_e)*dx/2;

            %___________________________________
            %unitary vectors normal to faces
            unorm_w=[-upar_w(2),upar_w(1)];
            unorm_n=[-upar_n(2),upar_n(1)];
            unorm_e=[upar_e(2),-upar_e(1)];
            unorm_s=[upar_s(2),-upar_s(1)];

            uVecNormFaces(i,j,:,:)=[unorm_w;unorm_n;unorm_e;unorm_s];

            %Vector from central node to center of faces
            %extract face centers coordinates
            %face_cent_w=reshape(faceCentrs(i,j,1,:),[1,2]);
            %face_cent_n=reshape(faceCentrs(i,j,2,:),[1,2]);
            %face_cent_e=reshape(faceCentrs(i,j,3,:),[1,2]);
            %face_cent_s=reshape(faceCentrs(i,j,4,:),[1,2]);

            %caculate vectors
            vec_p_fw=fcent_w-cell_centroid;
            vec_p_fn=fcent_n-cell_centroid;
            vec_p_fe=fcent_e-cell_centroid;
            vec_p_fs=fcent_s-cell_centroid;
            

            %Store vectors
            VecCentNodFaceCents(i,j,:,:)=[vec_p_fw;vec_p_fn;vec_p_fe;...
                vec_p_fs];


        end
    end

    %Variables that require previously computed values 


    %INTERIOR CELLS 
    %distance between central node and neigborhood nodes
    %unitary vectors from central node to ngb nodes
   

    for i=2:ny-1 %Iterate over rows
        for j=2:nx-1 %Iterate over columns

            %get the centroids of the nodes 
            cent_p=reshape(cellCentrs(i,j,:),[1,2]);
            cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
            cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
            cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
            cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

            vec_p_w=cent_w-cent_p;
            vec_p_n=cent_n-cent_p;
            vec_p_e=cent_e-cent_p;
            vec_p_s=cent_s-cent_p;

            %Get the distances between nodes 
            lgt_p_w=vecnorm(vec_p_w);
            lgt_p_n=vecnorm(vec_p_n);
            lgt_p_e=vecnorm(vec_p_e);
            lgt_p_s=vecnorm(vec_p_s);

            %STORE the  distances
            distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


            %Calculate unitary vectors 
            uvec_p_w=vec_p_w/lgt_p_w;
            uvec_p_n=vec_p_n/lgt_p_n;
            uvec_p_e=vec_p_e/lgt_p_e;
            uvec_p_s=vec_p_s/lgt_p_s;

            %STORE the values 

            uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

            %WEIGHTED LEAST SQUARE OPERATOR

            %Neighborhood nodes grouped in a vector
            ngb_nodes=[cent_w;cent_n;cent_e;cent_s];

            dist_pn=ngb_nodes - cent_p;
            %Vector of weights
            norm=vecnorm(dist_pn,2,2);
            w_v=1./norm;
            w_v=diag(w_v);
            %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
            G_m=dist_pn'*(w_v*dist_pn);
            gradient_weights_op=G_m\(dist_pn'*w_v);

            wlsqOperator(i,j,:,:)=gradient_weights_op;
         
        end
    end

   
    %If there is no neighborhood node is calculated
    %the distance to the center of the face instead

    %Sites where no neighboorhood nodes available
    
    %left  wall (west)
    j=1;
    for i =2:ny-1
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(faceCentrs(i,j,1,:),[1,2]);%face center
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];
       
    end

    %Top wall north
    i=1;
    
    for j=2:nx-1
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(faceCentrs(i,j,2,:),[1,2]);%face center
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    end

    %1.2.3- Right Wall - East 
    
    j=nx;
    for i=2:ny-1

        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(faceCentrs(i,j,3,:),[1,2]);%face center
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];        

    end

    %1.2.4- Bottom wall -- South
    
    i=ny;
    
    for j=2:nx-1
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(faceCentrs(i,j,4,:),[1,2]);%face center

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    end

    %CORNERS
    %1.3.1- North-West
    i=1;
    j=1;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(faceCentrs(i,j,1,:),[1,2]);%face center
    cent_n=reshape(faceCentrs(i,j,2,:),[1,2]);%face center
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %1.3.2- North-East
    i=1;
    j=nx;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=reshape(faceCentrs(i,j,2,:),[1,2]);%face center
    cent_e=reshape(faceCentrs(i,j,3,:),[1,2]);%face center
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %1.3.3- South-East (Symetric at South and outlet at East)
    i=ny;
    j=nx;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=reshape(faceCentrs(i,j,3,:),[1,2]);%face center
    cent_s=reshape(faceCentrs(i,j,4,:),[1,2]);%face center

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %1.3.4- South-West
    i=ny;
    j=1;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(faceCentrs(i,j,1,:),[1,2]);%face center
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=reshape(faceCentrs(i,j,4,:),[1,2]);%face center

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %______________________________________________________________
    %Additional points in faces for weighted least square gradients
    %West Edge 

    %aditional centers for least square gradient west edge
    centPointAddBoundWest=zeros(ny,2);
    %aditional centers for least square gradient North edge
    centPointAddBoundNorth=zeros(nx,2);
    %aditional centers for least square gradient East edge
    centPointAddBoundEast=zeros(ny,2);
    %aditional centers for least square gradient south edge
    centPointAddBoundSouth=zeros(nx,2);


    j=1;
    for i=1:ny
        
        centPointAddBoundWest(i,:)=reshape(faceCentrs(i,j,1,:),[1,2]);
    end
    %North
    i=1;
    for j=1:nx
        centPointAddBoundNorth(j,:)=reshape(faceCentrs(i,j,2,:),[1,2]);
    end
    %East
    j=nx;
    for i=1:ny
        centPointAddBoundEast(i,:)=reshape(faceCentrs(i,j,3,:),[1,2]);
    end
    %South (Special treatment to simulate a perpendicular line to the
    %solid surface
    i=ny;
    for j =1:nx

        %Get vertex west south from cell
        vert_ws= reshape(cellVertxs(i,j,4,:),[1,2]);
        %get cell center
        cel_cent=reshape(cellCentrs(i,j,:),[1,2]);
        %get unitary vector tangential to face S
        eta_s =reshape(uVecParlFaces(i,j,4,:),[1,2]);

        %calculate point coordinates for wall pressure boundary
        %the line that joints this point to cell center is perpendicular 
        %to face a, then least weighted squares gradient can be used to 
        % compute
        %pressure gradient at wall with non orthogonal faces
        centPointAddBoundSouth(j,:)=vert_ws +...
            dot((cel_cent-vert_ws),eta_s)*eta_s;

    end

    %Replace face center at south for centPointAddBoundSouth
    i=ny;
    for j=1:nx
        faceCentrs(i,j,4,:)=centPointAddBoundSouth(j,:);
    end

    %Additional Weighted Least Square Gradients

    %west edge
    j=1;
    for i = 2:ny-1
        %get the centroids of the nodes (replace west node to additional
        %center west face node
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=centPointAddBoundWest(i,:);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        %WEIGHTED LEAST SQUARE OPERATOR

        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];

        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);

        wlsqOperator(i,j,:,:)=gradient_weights_op;

    end
    
    %North Edge
    i=1;
    for j=2:nx-1
        %node center coordinates
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=centPointAddBoundNorth(j,:);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
        
        %WEIGHTED LEAST SQUARE OPERATOR
        
        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
        
        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);
        
        wlsqOperator(i,j,:,:)=gradient_weights_op;

    end
    
    %East Edge 
    j=nx;
    for i =2:ny-1
        %node center coordinates
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=centPointAddBoundEast(i,:);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
        
        %WEIGHTED LEAST SQUARE OPERATOR
        
        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
        
        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);
        
        wlsqOperator(i,j,:,:)=gradient_weights_op;


    end

    %South Edge 

    i=ny;
    for j=2:nx-1
        %node center coordinates
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=centPointAddBoundSouth(j,:);
        
        %WEIGHTED LEAST SQUARE OPERATOR
        
        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
        
        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);
        
        wlsqOperator(i,j,:,:)=gradient_weights_op;

    end

    %West north corner
    i=1;
    j=1;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=centPointAddBoundWest(i,:);
    cent_n=centPointAddBoundNorth(j,:);
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;



    %east north corner
    i=1;
    j=nx;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=centPointAddBoundNorth(j,:);
    cent_e=centPointAddBoundEast(i,:);
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;



    %east south corner
    i=ny;
    j=nx;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=centPointAddBoundEast(i,:);
    cent_s=centPointAddBoundSouth(j,:);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;

    %west south corner
    i=ny;
    j=1;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=centPointAddBoundWest(i,:);
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=centPointAddBoundSouth(j,:);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;

    %_______________________vertexes vector_______________________________

    %vectors from cell center each to vertex
     for i=1:ny
        for j=1:nx

            %get cell vertex
            cell_vertex=reshape(cellVertxs(i,j,:,:),[4,2]);
            vrt_wn=cell_vertex(1,:);
            vrt_en=cell_vertex(2,:);
            vrt_es=cell_vertex(3,:);
            vrt_ws=cell_vertex(4,:);

            %get cell centroid
            cell_cent=reshape(cellCentrs(i,j,:),[1,2]);

            %compute vectors;
            vec_c_wn=vrt_wn-cell_cent;
            vec_c_en=vrt_en-cell_cent;
            vec_c_es=vrt_es-cell_cent;
            vec_c_ws=vrt_ws-cell_cent;

            %Store values 

            VecCentVertx(i,j,:,:)=[vec_c_wn;vec_c_en;vec_c_es;vec_c_ws];

        end
     end

     %_________________Weight factors____________________________________

     for i=1:ny
         for j=1:nx

             %get distances from cell center to neighborhood nodes
             dist_nods=reshape(distNeighbNods(i,j,:,:),[1,4]);

             %split dist_nods
             dist_p_w=dist_nods(1);
             dist_p_n=dist_nods(2);
             dist_p_e=dist_nods(3);
             dist_p_s=dist_nods(4);

             %get distances from cell center to each face center

             %get vectors from cell center to face centers

             vec_cent_faces=reshape(VecCentNodFaceCents(i,j,:,:),[4,2]);

             %split vec_cent_faces

             vec_cent_w=vec_cent_faces(1,:);
             vec_cent_n=vec_cent_faces(2,:);
             vec_cent_e=vec_cent_faces(3,:);
             vec_cent_s=vec_cent_faces(4,:);

             %get the lenght of the vectors

             dist_cent_fw=vecnorm(vec_cent_w);
             dist_cent_fn=vecnorm(vec_cent_n);
             dist_cent_fe=vecnorm(vec_cent_e);
             dist_cent_fs=vecnorm(vec_cent_s);

             %compute weight factors
             weight_c_w=(dist_p_w - dist_cent_fw)/dist_p_w;
             weight_c_n=(dist_p_n - dist_cent_fn)/dist_p_n;
             weight_c_e=(dist_p_e - dist_cent_fe)/dist_p_e;
             weight_c_s=(dist_p_s - dist_cent_fs)/dist_p_s;

             %Store values

             weightDistFactors(i,j,:)=[weight_c_w,weight_c_n,weight_c_e...
                 ,weight_c_s];

         end
     end

     %___________________Mesh with cell centers__________________________
    for i=1:ny
        for j=1:nx
            node_cord=reshape(cellCentrs(i,j,:),[1,2]);
            Xctrs(i,j)=node_cord(1);
            Yctrs(i,j)=node_cord(2);

        end
    end
    %___________MINIMUM DISTANCE TO THE WALL COMPUTATION_____________________
   
    %k_1 = pi/0.9;
    %k_2 = k_1*(distPlat + 0.3);
   
    bump_start = distPlat + 0.3;
    bump_end = distPlat + 1.2;
    plate_end = distPlat + 1.5;  % Assuming bumpLgt = 1.5
    %{
    
    
    % Define the bump function and its derivatives as function handles
    bump_f = @(x) 0.05 * (sin(k_1*x - k_2)).^4;
    bump_df = @(x) 0.2 * k_1 * (sin(k_1*x - k_2)).^3 .* cos(k_1*x - k_2);
    bump_d2f = @(x) 0.6 * k_1^2 * (sin(k_1*x - k_2)).^2 .* ...
        (cos(k_1*x - k_2).^2 - sin(k_1*x - k_2).^2);
    
    % Define the objective function and its derivative for Newton-Raphson
    g = @(x, x_cent, y_cent) 2*(x - x_cent) + 2*(bump_f(x) - y_cent) .* bump_df(x);
    dg = @(x, x_cent, y_cent) 2 + 2*(bump_df(x).^2 + (bump_f(x) - y_cent) .* bump_d2f(x));
    %}
    for i = 1:ny
        for j = 1:nx
            % Get the x & y coordinate of the center of the cell
            x_cent = cellCentrs(i,j,1);
            y_cent = cellCentrs(i,j,2);
            
            % Case 1: Point above the bump region
            if x_cent > bump_start && x_cent < bump_end
                %if i > ny-12


                    %cells very close to the wall

                    %use dot product aproximation

                    %get south west vertex of the cell ny,j

                    vert_sw=reshape(cellVertxs(ny,j,4,:),[1,2]);

                    %get south face normal unitary vector from cell ny,j

                    norm_s= reshape(uVecNormFaces(ny,j,4,:),[1,2]);

                    %vector from vertex sw to cell centroid

                    vec_cent_sw =[x_cent,y_cent] - vert_sw;

                    %dot vector with inverse unitary normal 

                    distMinWall(i,j)= dot(vec_cent_sw,(-norm_s)); 
                %{
                else
                    % Use Newton-Raphson to find the closest point on the bump
                    
                    % Improved initial guess - use the x-coordinate of the centroid
                    % but clamp it to the bump region
                    x_init = max(bump_start, min(x_cent, bump_end));
                    
                    % Newton-Raphson with convergence criteria
                    max_nr_iter = 100;  % Reduced from 1000
                    tolerance = 1e-12;
                    
                    for m = 1:max_nr_iter
                        g_val = g(x_init, x_cent, y_cent);
                        dg_val = dg(x_init, x_cent, y_cent);
                        
                        % Avoid division by zero and check convergence
                        if abs(dg_val) < 1e-12
                            % Derivative too small, use alternative method
                            x_init = (bump_start + bump_end) / 2;
                            break;
                        end
                        
                        x_new = x_init - g_val / dg_val;
                        
                        % Clamp to bump region to stay within valid domain
                        x_new = max(bump_start, min(x_new, bump_end));
                        
                        % Check convergence
                        if abs(x_new - x_init) < tolerance
                            x_init = x_new;
                            break;
                        end
                        
                        x_init = x_new;
                        
                        % Additional safeguard: if we're not making progress
                        if m > 10 && abs(g_val) > 1e9
                            % Fallback: use centroid's x projected to bump region
                            x_init = max(bump_start, min(x_cent, bump_end));
                            break;
                        end
                    end
                    
                    % Calculate the minimum distance
                    f_x = bump_f(x_init);
                    distMinWall(i,j) = sqrt((x_cent - x_init)^2 + (y_cent - f_x)^2);
                end
                %}
            % Case 2: Point above the flat portion of the plate
            elseif (x_cent >= distPlat && x_cent < bump_start) || ...
                   (x_cent >= bump_end && x_cent < plate_end)
                
                distMinWall(i,j) = y_cent;
                    
            % Case 3: Cell centroid before or after the plate
            else
                if x_cent < distPlat
                    % Point before the plate - distance to leading edge
                    p_wall = [distPlat, 0];
                else
                    % Point after the plate - distance to trailing edge
                    p_wall = [plate_end, 0];
                end
                distMinWall(i,j) = sqrt((x_cent - p_wall(1))^2 + (y_cent - p_wall(2))^2);
            end
        end
    end


end