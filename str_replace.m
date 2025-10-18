%FOR SEARCH AND REPLACE AND BRIEF CODES TEST

%SKEWNESS TEST 
%{
dot_test=zeros(1,nx);

i=40;
for j=1:nx
    %get vector normal to face e
    normal_fe=reshape(uVecNormFaces(i,j,3,:),[1,2]);
    %get vector from node p to node E
    vector_p_e = reshape(uVecNeighbNods(i,j,3,:),[1,2]);
    % Calculate the dot product between the normal vector and the vector from p to E
    dot_test(j) = dot(normal_fe, vector_p_e);

    % Store the skewness angle based on the dot product
    %skewness_angle(j) = acos(dot_test(j));


end

% Store the skewness angle based on the dot product wtih degrees
skewness_angle=acosd(dot_test);
%}

%Get  all vertex for cell [31,38]
%{
i=31; j=38;

cel_vert = reshape(cellVertxs(i,j,:,:),[4,2]);
vert_wn=cel_vert(1,:);
vert_en=cel_vert(2,:);
vert_es=cel_vert(3,:);
vert_ws=cel_vert(4,:);

%get cell center
cel_cent=reshape(cellCentrs(i,j,:),[1,2]);

%get unitary vector tangential to face S
eta_s =reshape(uVecParlFaces(i,j,4,:),[1,2]);

%calculate point coordinates for wall pressure boundary
%the line that joints this point to cell center is perpendicular 
%to face a, then least weighted squares gradient can be used to compute
%pressure gradient at wall with non orthogonal faces

cent_ps_wallBC= vert_ws +  dot((cel_cent-vert_ws),eta_s)*eta_s;


%plot as points

figure(2)
hold on
%a red circle ('ro')
plot(vert_wn(1), vert_wn(2), 'ro', 'MarkerSize', 8,'MarkerFaceColor', 'r');
plot(vert_en(1), vert_en(2), 'ro', 'MarkerSize', 8,'MarkerFaceColor', 'r');
plot(vert_es(1), vert_es(2), 'ro', 'MarkerSize', 8,'MarkerFaceColor', 'r');
plot(vert_ws(1), vert_ws(2), 'ro', 'MarkerSize', 8,'MarkerFaceColor', 'b');

% Plot cell center
plot(cel_cent(1), cel_cent(2), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
%plot artifitial point to simulate presure gradient zero perpendicular to
%wall
plot(cent_ps_wallBC(1), cent_ps_wallBC(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

hold off
grid on

%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

%node center coordinates
cent_p=reshape(cellCentrs(i,j,:),[1,2]);
cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
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

%aditional centers for least square gradient west edge
centPointAddBoundWest=zeros(ny,2);
%aditional centers for least square gradient North edge
centPointAddBoundNorth=zeros(nx,2);
%aditional centers for least square gradient East edge
centPointAddBoundEast=zeros(ny,2);
%aditional centers for least square gradient south edge
centPointAddBoundSouth=zeros(nx,2);


%___________________________________________
orthogonal_Test=zeros(1,nx);
i=ny;
for j=1:nx

    %get face center
    center_test= reshape(cellCentrs(i,j,:),[1,2]);

    %get west south vertex
    vws=reshape(cellVertxs(i,j,4,:),[1,2]);

    %get face center at fase s
    fs_cent=reshape(faceCentrs(i,j,4,:),[1,2]);


    %vector from vertex ws to face s center
    v_1=fs_cent-vws;

    %vector from face s center to  cell center

    v_2=center_test-fs_cent;

    % Calculate the dot product to check orthogonality
    orthogonal_Test(j) = acosd(dot(v_1, v_2));



end
%}
%Skewness test for north faces 

%create a matrix to store values

%{

skw_northFc=zeros(ny,nx);

for i=1:ny
    for j=1:nx
        %Get vector from cell center to north node
        vec_c_n = reshape(uVecNeighbNods(i,j,2,:),[1,2]);

        %get vector normal to face N

        vec_fn=reshape(uVecNormFaces(i,j,2,:),[1,2]);

        skw_northFc(i,j) = acosd(dot(vec_c_n,vec_fn));

    end
end
%}

%correction of pressure correction code

%(fAw/dist_w)*((volP/aP(i,j))*weight_w  + (1-weight_w)*(volW/aP(i,j-1)))

%pressure coeffitient code 

%x_pressure_coeff=zeros(size(nxSolid));

%for j=1:nx_fp
%    j_1=nxSolid(j);
%    x_pressure_coeff(j)=reshape(faceCentrs(ny,j_1,4,1),[1,1]);
%end 
%x_pressure_coeff=x_pressure_coeff-(domLgt-bumpLgt)/2;
%pressure_coeff=(p(ny,nxSolid)-p0)/rho*u0;

%figure(13)
%plot(x_pressure_coeff,pressure_coeff ,'- o')
%title("PRESSURE COEFFITIENT ON THE PLATE ")
%xlabel("X POSITION")
%ylabel("PRESSURE COEFFITIENT ")

%Vorticity
%vorticity=zeros(ny,nx);
%{
for i=1:ny
    for j=1:nx
        dv_dx=reshape(grad_v(i,j,1,1),[1,1]);
        du_dy=reshape(grad_u(i,j,1,2),[1,1]);
        vorticity(i,j)=dv_dx-du_dy;           
    end 
end

figure(8)
contourf(Xctrs, Yctrs, vorticity, 20, 'LineColor', 'none');
colorbar
title('Vorticity Field')
xlabel('X')
ylabel('Y')
axis equal
colormap jet
%}


%minimum distance to the walll

%{
for i=1:nx
    for j= 1:ny

        %Get the x & y coordinate of the center of the cell
        x_cent=reshape(cellCentrs(i,j,1),[1,1]);
        y_cent=reshape(cellCentrs(i,j,2),[1,1]);

        %case 1) Node point lies backwars or upwards plate
        if x_cent < distPlat || x_cent > distPlat + 1.5
            if x< distPlat
                p_wall=[distPlat + 0.7500,0];
                d_wall=sqrt(( x - p_wall(1) )^2 + (y - p_wall(2))^2 );
            else
                p_wall=[distPlat + 1.5,0];
                d_wall=sqrt(( x - p_wall(1) )^2 + (y - p_wall(2))^2 );       
            end
        elseif  x < distPlat + 0.3 || x > distPlat + 1.2
            d_wall=y;    
        else
        %case 2) Node is above the plate    
        %use  Newton-Rhapson to find nearest distance to the wall
        %find nearest point using newton raphson
        [x_n,y_n] = Newton_Raphson_nd(x,y);
        %case 2.1) calculate with newton Raphson
        d_walln=sqrt(( x - x_n)^2 + (y - y_n)^2 );
        %case 2.2) just take the y coordinate as the distance 
        d_wall_b= y;
        d_wall=min(d_wall_b,d_walln);%compare 
        end

    end
end
%}
%___________MINIMUN DISTANCE TO THE WALL COMPUTATION_____________________

%distPlat= distance from origin to the leading edge of the plate
%k_1=pi/0.9;% constant 1 for defining the bump 

%k_2=k_1*(distPlat + 0.3);% constant 2  for defining the bump 


%Bump is described by the formula f(x) = 0.05(sin(k_1*x - k_2)^4
%from distPlat + 0.3 < x< distPlat + 1.2 
%square distance from any point (x_cent,y_cent) to a point in the bump is
%d^2 = (x - x_cent)^2 + (f(x) - y_cent)^2
%Computing the derivate and equaling to zero to find the minimum distance
%0=2*(x- x_cent) + 2*(f(x) - y_cent)*f'(x) =g(x)
%g'(x)=2 + 2*((f'(x))^2 + (f(x) -y_cent)*f'(x))
%Newton Rhapson formula 
%xn+1=x_n + g(x_n)/g'(x_n)
%{
for i=1:ny
    for j=1:nx

        %Get the x & y coordinate of the center of the cell
        x_cent=reshape(cellCentrs(i,j,1),[1,1]);
        y_cent=reshape(cellCentrs(i,j,2),[1,1]);

        %case 1: Point above the bump
        if x_cent > (distPlat + 0.3 ) && x_cent < (distPlat + 1.2)
            %Use Newton Rhapson
            %Newton - Raphson method to obtain calculate nearest distance to the wall
           
            %Initial Guess
            x_init=distPlat + 0.3;
            
            for m=1:1000
                f_x=0.05*((sin(k_1*x_init -k_2))^4);
                df_dx=0.2*k_1*((sin(k_1*x_init -k_2))^3)*cos(k_1*x_init -k_2);%f'(x)
                df2_dx2=0.6*(k_1^2)*(((sin(k_1*x_init -k_2))^2)*cos(k_1*x_init -k_2) - ((sin(k_1*x_init -k_2))^4));%f''(x)
                g_x=2*(x_init -x_cent) + 2*(f_x - y_cent)*df_dx;%g(x)
                d2g_dx2= 2 + 2*(df_dx^2 + (f_x -y_cent)*df2_dx2);%g´(x)
                %_______MAIN NEWTON RAPHSON_______________________
                x_init =x_init - g_x/d2g_dx2;
            
            end
                        
            distMinWall(i,j)=vecnorm([x_cent,y_cent]-[x_init,f_x]);

        %case 2: point above the flat zone of the plate
        elseif x_cent>=distPlat && x_cent<(distPlat + 0.3 ) || x_cent > (distPlat + 1.2) && x_cent < distPlat + bumpLgt

            %minimum distance is equal to the 
            distMinWall(i,j)=y_cent;
                
        %case 1: Cell centroid backward or upward the plate
        else
            %point backward plate
            if x_cent < distPlat
                p_wall=[distPlat,0];
                distMinWall(i,j)=vecnorm([x_cent,y_cent]-p_wall);

            %point upward plate     
            else
                p_wall=[distPlat+bumpLgt,0];
                distMinWall(i,j)=vecnorm([x_cent,y_cent]-p_wall);
                
            end 
        end

    end
end
%}

%figure(11)
%contourf(Xctrs(ny:-1:ny-20,nxBump),Yctrs(ny:-1:ny-20,nxBump),distMinWall(ny:-1:ny-20,nxBump))
%title("Minimum Wall Distance")
%xlabel("Height (m)")
%ylabel("Lenght (m)")
%colormap jet
%colorbar
%axis equal


%sqrinversedistMinWall=(1./distMinWall.^2);

%figure(12)
%contourf(Xctrs(ny:-1:ny-27,nxBump),Yctrs(ny:-1:ny-27,nxBump),sqrinversedistMinWall(ny:-1:ny-27,nxBump))
%title("Minimum inverse squared Wall Distance")
%xlabel("Height (m)")
%ylabel("Lenght (m)")

%{
%}

%{
figure(11)
contourf(Xctrs,Yctrs,distMinWall)
title("Minimum Wall Distance")
xlabel("Height (m)")
ylabel("Lenght (m)")
colormap jet
colorbar
axis equal


sqrinversedistMinWall=(1./distMinWall.^2);

figure(12)
contourf(Xctrs,Yctrs,sqrinversedistMinWall)
title("Minimum Wall Distance")
xlabel("Height (m)")
ylabel("Lenght (m)")


%

%-cellVols(i,j)*(grad_px + grad_tau_xx_x + grad_tau_xy_y)
%-cellVols(i,j)*(grad_py + grad_tau_xy_x + grad_tau_yy_y)
%}

%evaluating y_cent 

%LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
%Implicit Taylor 1st order expansion
%Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
% - nu_tilde_k)
%central coeffitient
%AP=-Q'(nu_tilde_k)*cellVol
%Source
%AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

%Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)
%d_wall=distMinWall(i,j);


%__________ Turbulence Production Tp(nu_tilde) ___________
i=27;
j=45;

d_wall=distMinWall(i,j);
%Viscous ratio
x_nu=nu_tilde_k/nu;
x_nu_prime=1/nu;

%%Viscous damping function 1 fv1
fv1=(x_nu^3)/(x_nu^3 + cv1^3);
fv1_prime=(x_nu_prime*3*(x_nu^2)*(cv1^3))/((x_nu^3 + cv1^3)^2);


%Viscous damping function 2 fv2
fv2=1 - x_nu/(1 + x_nu*fv1);
fv2_prime=(x_nu^2*fv1_prime - x_nu_prime)/(1 + x_nu*fv1)^2;
 
%abs  vorticity (gamma tilde)
s_vort=abs(vorticity(i,j));

%modified vorticity

s_vort_tilde = s_vort + (nu_tilde_k*fv2/(kappa^2 * d_wall^2));
s_vort_tilde_prime= (1/(kappa^2 * d_wall^2))*(fv2 + fv2_prime*nu_tilde_k);

%tp Turbulence production
tp=cb1*s_vort_tilde*nu_tilde_k;
tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);

%__________ Turbulence Destruction Td(nu_tilde) __________

%function r
r=min(10,nu_tilde_k/(s_vort_tilde * kappa^2 * d_wall^2 ));
r_prime=(1/((s_vort_tilde^2)*(kappa^2 * d_wall^2)))*( s_vort_tilde -nu_tilde_k*s_vort_tilde_prime);

%function g
g=r + cw2*(r^6 -r);
g_prime=r_prime*(1-cw2 + 6*cw2*r^5);

%wall damping function fw
cw3_6=cw3^6;


fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
fw_prime=g_prime*(cw3_6)*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

%Td
td=-cw1*fw*(nu_tilde_k/d_wall)^2;
td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + 2*fw*(nu_tilde_k/d_wall));

figure(12)
contourf(Xctrs(ny:-1:ny-27,nxBump),Yctrs(ny:-1:ny-27,nxBump),mu_turbulent(ny:-1:ny-27,nxBump))
title("Minimum inverse squared Wall Distance")
xlabel("Height (m)")
ylabel("Lenght (m)")