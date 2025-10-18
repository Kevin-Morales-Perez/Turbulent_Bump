function [tau_xx,tau_xy,tau_yy] = saBoussinesqTurbulentStresses(...
    mu_turbulent,grad_u,grad_v,tau_xx,tau_xy,tau_yy,nx,ny)
    %COMPUTATION OF EDDY VISCOSITY FROM SPALART - ALLMARAS MODEL 
    % TAU_ij=MU_T*2*(S_ij)

    for i=1:ny
        for j=1:nx
            
            %Split gradients in componets

            du_dx=grad_u(i,j,1,1);  
            du_dy=grad_u(i,j,1,2);

            dv_dx=grad_v(i,j,1,1);
            dv_dy=grad_v(i,j,1,2);

            tau_xx(i,j)= mu_turbulent(i,j)*(2*du_dx);
            tau_xy(i,j)= mu_turbulent(i,j)*(du_dy + dv_dx);
            tau_yy(i,j)= mu_turbulent(i,j)*(2*dv_dy);
        end
    end
end