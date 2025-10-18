function [mu_turbulent] = saEddyViscosity(nu_tilde,nu,rho,cv1)
    %COMPUTATION OF EDDY VISCOSITY ACCORDING TO SPALART - ALLMARAS MODEL
    %WITH DAMPING FUNCTIONS
   
        X_sa=nu_tilde./nu;%intermediate variable

        fv1= (X_sa.^3)./(X_sa.^3 + cv1^3); %Empirical viscous damping function in turbulence modell
        
        nu_turbulent=max(0,nu_tilde.*fv1);%Kinematic eddy viscosity

        mu_turbulent=rho.*nu_turbulent;%Molecular Eddy viscosity
    
end

