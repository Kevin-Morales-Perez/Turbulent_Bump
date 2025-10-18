%BUMP  MESH ORIGINAL
%CODE  TO GENERATE MESH OVER A BUMP WITH INFLATION LAYERS AT THE 
%BOUNDARY, LEADING AND TRAILING EDGES 
%THERE IS A FREE - STREAM ZONE AT THE LEADING EDGE AND OTHER IN THE 
%TRAILING EDGE 


close all
%% DATA

%Geometry
bumpLgt=1.5;                         %bump length (m)
domHgt=3;                            %Domain Height (m)
domLgt=3;                          %Domain Length (m)

distPlat=0.5*(domLgt - bumpLgt);     %Distance from origin to leading edge
%Should be always placed at the center 

if (distPlat+ bumpLgt)<domLgt
    %VALID GEOMETRY
    
    %NUMBER OF CELLS IN EACH SECTION
    nx=90;% number of cells in x axis 
    ny=30;% number of cells in y axis
    ncells=nx*ny; %Total cells

    x=zeros(1,nx+1); %grid positions for x
    y=zeros(1,ny+1); %grid positions for y

    %Calculate the ratio between the length of the plate and the
    %  domain length
    fpRatio = 0.5;%bumpLgt/domLgt;% Now the plate should contain half o f the 
    %total cells in X

    %Calculate the ratio between the distance to the leading edge and the
    %  domain lenght origin (we will call this upstream zone)
    
    upstrRatio=0.5*(1-fpRatio);%distPlat/domLgt; %%for 50%


    %Divide proportionaly the number of cells in x between the plate 
    % and the remaining space 

    nx_fp=round(nx*fpRatio,0);%cell in the plate 
    nx_upstr=round(nx*upstrRatio);%Cells in the upstream zone 
    nx_dwnstr=nx-(nx_fp + nx_upstr);%Cells in the downstream zone

    nxSymcond=[1:1:nx_upstr,nx-nx_dwnstr+1:1:nx];%indexes for symetric 
    nxSolid=nx_upstr+1:nx-nx_dwnstr; %Indexes for solid plate 
    %boundary condition

    %INFLATION LAYERS

    %INFLATION LAYER FOR X AXIS 

    dx_1=0.00331*3;%Lenght of the first cell at the leading and trailing edge
    grthFx=1.3;% Geometrical growth  Factor
    ninfLx=2;%Cells in the l/t edges side of the inflation layers

    %Calculate grid spaces for x in the inflation layer

    %define where  finishes the upstream zone and begins the plate
    x(nx_upstr +1)=distPlat;

    %Define where finishes the flat plate
    x(nx_upstr+nx_fp+1)=distPlat+ bumpLgt;

    %Define where finishes the domain
    x(nx+1)=domLgt;

    %Calculate grid spaces for x in the inflation layer (leading edge ,RS)
    %At the same time mirror the steps for the free stream portion
    %repeat for the trailing edge
    for i=1:ninfLx
        dx_inf=dx_1*grthFx^(i-1);%Layer lenght
        x(nx_upstr+1+i)=dx_inf + x(nx_upstr+i);%Flat plate LE
        x(nx_upstr+1-i)=-dx_inf +x(nx_upstr+2-i);%upwind portion
        
        %Apply same for the trailing edge

        x(nx_upstr + nx_fp + 1 + i) = dx_inf + ...
            x(nx_upstr + nx_fp + i);%Flat plate TE
        x(nx_upstr + nx_fp + 1-i) = -dx_inf + ...
            x(nx_upstr + nx_fp + 2-i);%Downstream portion
    end

    %SMOTH TRANSITIONS FOR X

    %SMOTH TRANSITION IN THE UPSTREAM ZONE

    dxMaxinf=x(nx_upstr + ninfLx +1)- x(nx_upstr + ninfLx); % Higher x grid space 
    % in inflation layers for X

    ninfLxupstr= nx_upstr - ninfLx;%Cells in the smooth transition 
    % upstream layer

    dxT2upstr=x(ninfLxupstr+1);%Lenght of the 
    %upstream smooth transition layer

    %%% Computig geometrical growth factor for upstream smooth 
    % transition layer

    grthFT2upstr=0.00001;%Geometrical grow factor 2
    dxT2upstrit=0;%Boundary layer height iterated
    while dxT2upstrit<=dxT2upstr
        grthFT2upstr=grthFT2upstr+0.00001;
        dxT2upstrit=dxMaxinf*((1-grthFT2upstr^ninfLxupstr)/...
            (1-grthFT2upstr));
    end

    %Fill the x grid vector with smoth transition downstream layer steps 
    
    for i=1:ninfLxupstr-1
        dx_inf=dxMaxinf*grthFT2upstr^(i-1);%layer lenght
        x(nx_upstr-ninfLx+1-i) = x(nx_upstr-ninfLx+2-i) - dx_inf;
    end


    %SMOTH TRANSITION IN THE DOWNSTREAM ZONE 

    ninfLxdwnstr=nx_dwnstr-ninfLx;%Cells in the smooth transition 
    % downstream layer

    dxT2dwnstr=domLgt-x(nx_upstr+nx_fp+ninfLx+1);%Lenght of the 
    %downstream smooth transition layer

    %%% Computig geometrical growth factor for downstream smooth 
    % transition layer

    grthFT2dwnstr=0.00001;%Geometrical grow factor 2
    dxT2dwnstrit=0;%Boundary layer height iterated
    while dxT2dwnstrit<=dxT2dwnstr
        grthFT2dwnstr=grthFT2dwnstr+0.00001;
        dxT2dwnstrit=dxMaxinf*((1-grthFT2dwnstr^ninfLxdwnstr)/...
            (1-grthFT2dwnstr));
    end

    %Fill the x grid vector with smoth transition downstream layer steps 
    
    for i=1:ninfLxdwnstr
        dx_inf=dxMaxinf*grthFT2dwnstr^(i-1);%layer lenght
        x(nx_upstr+nx_fp+ninfLx+1+i) = x(nx_upstr+nx_fp+ninfLx+i) + dx_inf;
    end


    %SMOTH TRANSITION IN THE PLATE

    %Calculate the distance from the origin to the midpoint of the plate
    distMidPlat=distPlat + 0.5*bumpLgt;

    %Smoth transition from leading edge inflation last layer to mid point
    %of the plate;

    %cells in the smooth transition zone of the plate (mid region)
    ninfLxmidFp=nx_fp-2*ninfLx;

    %Divide this distance in a upstream and downstream zone
    ninfLxmidFupstr= round(0.5*ninfLxmidFp);% Cells in the upstream zone

    ninfLxmidFdwnstr=ninfLxmidFp-ninfLxmidFupstr;%Cells in the downstream z.

    %SMOTH TRANSITION IN THE UPSTREAM ZONE (PLATE)

    dxT2upstrmidPlat=distMidPlat - x(nx_upstr + ninfLx + 1);%Lenght of the 
    %upstream smooth transition layer

    %%% Computig geometrical growth factor for upstream smooth 
    % transition layer

    grthFT2upstrmidPlat=0.00001;%Geometrical grow factor 2
    dxT2upstritmidPlat=0;%Boundary layer height iterated
    while dxT2upstritmidPlat<=dxT2upstrmidPlat
        grthFT2upstrmidPlat=grthFT2upstrmidPlat+0.00001;
        dxT2upstritmidPlat=dxMaxinf*((1-grthFT2upstrmidPlat^ninfLxmidFupstr)/...
            (1-grthFT2upstrmidPlat));
    end

    %Fill the x grid vector with smoth transition upstream layer steps 
    %for the plate 
    
    for i=1:ninfLxmidFupstr
        dx_inf=dxMaxinf*grthFT2upstrmidPlat^(i-1);%layer lenght
        x(nx_upstr + ninfLx + 1 + i) = x(nx_upstr + ninfLx + i)+ dx_inf;
    end

    %SMOTH TRANSITION IN THE DOWNSTREAM ZONE (PLATE)

    %Smoth transition from trailing edge inflation last layer to mid point
    %of the plate;

    %Lenght of the downstream smooth transition layer

    dxT2dwnstrmidPlat=x(nx - (nx_dwnstr + ninfLx) +1) - distMidPlat;

    %%% Computig geometrical growth factor for downstream smooth 
    % transition layer on plate

    grthFT2dwnstrmidPlat=0.00001;%Geometrical grow factor 2
    dxT2dwnstritmidPlat=0;%Boundary layer height iterated
    while dxT2dwnstritmidPlat<=dxT2dwnstrmidPlat
        grthFT2dwnstrmidPlat=grthFT2dwnstrmidPlat+0.00001;
        dxT2dwnstritmidPlat=dxMaxinf*((1-grthFT2dwnstrmidPlat^ninfLxmidFdwnstr)/...
            (1-grthFT2dwnstrmidPlat));
    end

    %Fill the x grid vector with smoth transition downstream layer steps 
    %for the plate 
    
    for i=1:ninfLxmidFdwnstr-1
        dx_inf=dxMaxinf*grthFT2dwnstrmidPlat^(i-1);%layer lenght
        x( nx - (nx_dwnstr + ninfLx) +1 - i) = x(nx - ...
            (nx_dwnstr + ninfLx) +2 - i)- dx_inf;
    end

    %dxVec=diff(x);%X grid spaces vector

    %Matrix for grid 
    X=zeros(ny+1,nx+1);
    for j=1:nx+1
        X(:,j)=x(j);
    end
    %______________   DEFINE THE BUMP POSITION __________________________

    %The original test case from the turbulence resource modeling stablish
    %the leading edge of the plate at the origin and the bump from 0.3 to 
    %1.2 then the upstream boundary is located in the fourth cuadrant of 
    %the cartesian plane, however here the upstream boundary begins at y=0
    %in such a way all the domain is in the first cuadrant 

    %find at which distance should begin the bump (or upstream)
    distupstrBump=distPlat+0.3;
    
    %find at which distance should finish the bump (or downstream)

    distdwnstrBump=distPlat+1.2;

    %find in the X grid vector the indexes at where the bump should start 
    % and finish

    bumpRangeIndices = find(x >= distupstrBump & x <= distdwnstrBump);

    %Extract from vector X the elements that correspont to the bump

    xBump=x(bumpRangeIndices);

    %calculate the y positions of the bump using the previos x coordinates

    x0Bump=xBump(1);%Take first x coordinate of the bump
    kBump=0.9/x0Bump; %Compute constant for the bump

    yBump=0.05*sin(pi*xBump/0.9 -(pi/kBump)).^4; %Y coordinates for bump

     %______________   Inflation layer parameters for Y axis  _____________

    dy_1=0.00167*0.25;% Height of the first cell (m);0.00021,0.0001;%
    grthF=1.35;% Geometrical growth  Factor ,1.3
    ninfL=10;% Cells in the boundary inflation layer

    %create the first row of the matrix Y , this row represents the 
    %shape and the boundary of the bump

    yFirstRow=zeros(1,nx+1);

    yFirstRow(bumpRangeIndices)=yBump;

    plot(x,yFirstRow);


    %create the matrix that will contain the grid of Y
    % at the moment as a test

    Ytest=zeros(ny+1,nx+1);
    
    %Replace the last row with the shape y vector for the bump
    Ytest(end,bumpRangeIndices)=yBump;

    ninfL2=ny-ninfL;%Cells in the second inflation layer for smooth transition

    for j=1:nx+1 %iterate over columns 

        %Create the boundary inflation layer
        for i=ny+1:-1:(ny+1-(ninfL-1))
            dy_inf=dy_1*grthF^(ny+1 -i);%layer lenght
            Ytest(i-1,j)=Ytest(i,j)+dy_inf;
        end

        %Calculating inflation layer parameters to fill the remaining
        %spaces in  y grid vector and have a smooth transition

        dy_2=Ytest(ny+1-ninfL,j)-Ytest(ny+2-ninfL,j); % Higher y grid space in boundary inf layer
        dy_2=dy_2*1.75; %increment max dy in the inf layer if needed
        dyt2=domHgt - Ytest(ny+1-ninfL,j);%Height of the second inflation layer

        %%% Computig geometrical growth factor for second inflation layer
        grthF2=0.00001;%Geometrical grow factor 2
        dyt2it=0;%Boundary layer height iterated
        while dyt2it<=dyt2
            grthF2=grthF2+0.00001;
            dyt2it=dy_2*((1-grthF2^ninfL2)/(1-grthF2));
        end

        % Fill the remaining y grid vector with second inflation layer
        for i=ny-ninfL:-1:2
            Ytest(i,j) = Ytest(i+1,j) + dy_2*grthF2^(ny-ninfL-i);
        end


    end

    Ytest(1,:)=domHgt;

    %Define in which x indexes is located the bump
    nxBump=[];
    for j=2:nx+1
        %SIMETRIC BC INDEXES 
        if Ytest(ny+1,j)>0
            nxBump=[nxBump,j-1];
        end

    end

    nxBump=[nxBump,(nxBump(end)+1)];


    %
    %Visualization of the grid.
    figure(1);
    % Plot vertical lines (constant j)
    plot(X, Ytest, 'Color', 'b', 'LineWidth', 0.5);
    hold on;
    % Plot horizontal lines (constant i)
    plot(X', Ytest', 'Color', 'b', 'LineWidth', 0.5);
    axis equal;
    xlabel('X');
    ylabel('Y');
    title(['Mesh for a flat plate with inflation layers at the' ...
        '     boundary and leading edge']);
    %}

    Y=Ytest;
    clear Ytest;

    dxVec=diff(x);

    %Save matrices X and Y 

    %save('bumpMesh.mat', 'X', 'Y');
   
else
    fprintf("DEFINE CORRECT VALUES FOR DOMAIN AND PLATE SIZES!!!")

end