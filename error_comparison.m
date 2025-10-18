%ERROR COMPARE CODE 

%1.- MOMENTUM LINK COEFFITIENTS AND SOURCES

aW_correct = aW;
aE_correct = aE;
aN_correct = aN;
aS_correct = aS;
aP_correct = aP;
aPv_correct = aPv;
suX_correct = suX;
suY_correct = suY;

%2.- SOLVE X MOMENTUM

% Store correct values
u_correct = u;
rsid_x_correct = rsid_x;
err_x_correct = err_x;


%3.- SOLVE Y MOMENTUM

% Store correct values
v_correct = v;
rsid_y_correct = rsid_y;
err_y_correct = err_y;

%4.- FACE VELOCITY COMPUTATION USING RIE - CHOW INTERPOLATION*

% Store correct values
u_face_correct = u_face;
v_face_correct = v_face;

%5.- PRESSURE CORRECTION LINK COEFFITIENTS AND MASS INBALANCE (SOURCE)

% Store correct values
ap_W_correct = ap_W;
ap_N_correct = ap_N;
ap_E_correct = ap_E;
ap_S_correct = ap_S;
ap_P_correct = ap_P;
suP_correct = suP;

%6.- SOLVE PRESSURE CORRECTION
% Store correct values
p_prime_correct = p_prime;
rsid_p_correct = rsid_p;
err_p_correct = err_p;


%7.- CORRECT PRESSURE

p_star_correct = p_star;

%8.- CORRECT CELL CENTER VELOCITY

% Store correct values
u_star_correct = u_star;
v_star_correct = v_star;

%9.- CORRECT FACE VELOCITY

% Store correct values
%u_face_correct = u_face;
%v_face_correct = v_face;


%10.- P=P_STAR
%p_correct = p;

%_______________________________  2  ____________________________________
%{
%1.- MOMENTUM LINK COEFFITIENTS AND SOURCES
% Compute errors for this step
error_aW = matrix_err_percentage(aW, aW_correct);
error_aE = matrix_err_percentage(aE, aE_correct);
error_aN = matrix_err_percentage(aN, aN_correct);
error_aS = matrix_err_percentage(aS, aS_correct);
error_aP = matrix_err_percentage(aP, aP_correct);
error_aPv = matrix_err_percentage(aPv, aPv_correct);
error_suX = matrix_err_percentage(suX, suX_correct);
error_suY = matrix_err_percentage(suY, suY_correct);

figure(1) ; surf(error_aW); title('aW Error'); colorbar;
figure(2) ;surf(error_aE); title('aE Error'); colorbar;
figure(3) ;surf(error_aN); title('aN Error'); colorbar;
figure(4) ;surf(error_aS); title('aS Error'); colorbar;
figure(5) ;surf(error_aP); title('aP Error'); colorbar;
figure(6) ;plot(error_aPv); title('aPv Error'); colorbar;
figure(7) ;surf(error_suX); title('suX Error'); colorbar;
figure(8) ;surf(error_suY); title('suY Error'); colorbar;

%2.- SOLVE X MOMENTUM

errors_u = matrix_err_percentage(u, u_correct);
errors_rsid_x = matrix_err_percentage(rsid_x, rsid_x_correct);
errors_err_x = matrix_err_percentage(err_x, err_x_correct);

figure(9); surf(errors_u); title('u Error'); colorbar;
figure(10); surf(errors_rsid_x); title('rsid_x Error'); colorbar;
figure(11); plot(errors_err_x); title('err_x Error'); 


%3.- SOLVE Y MOMENTUM

% Compute errors
errors_v = matrix_err_percentage(v, v_correct);
errors_rsid_y = matrix_err_percentage(rsid_y, rsid_y_correct);
errors_err_y = matrix_err_percentage(err_y, err_y_correct);

% Plot errors

figure(12); surf(errors_v); title('v Error'); colorbar;
figure(13); surf(errors_rsid_y); title('rsid_y Error'); colorbar;
figure(14); plot(errors_err_y); title('err_y Error');

%4.- FACE VELOCITY COMPUTATION USING RIE - CHOW INTERPOLATION*

% Compute errors
errors_u_face = matrix_err_percentage(u_face, u_face_correct);
errors_v_face = matrix_err_percentage(v_face, v_face_correct);

% Plot errors

figure(15); surf(errors_u_face); title('u\_face Error'); colorbar;
figure(16); surf(errors_v_face); title('v\_face Error'); colorbar;

%5.- PRESSURE CORRECTION LINK COEFFITIENTS AND MASS INBALANCE (SOURCE)

% Compute errors
errors_ap_W = matrix_err_percentage(ap_W, ap_W_correct);
errors_ap_N = matrix_err_percentage(ap_N, ap_N_correct);
errors_ap_E = matrix_err_percentage(ap_E, ap_E_correct);
errors_ap_S = matrix_err_percentage(ap_S, ap_S_correct);
errors_ap_P = matrix_err_percentage(ap_P, ap_P_correct);
errors_suP = matrix_err_percentage(suP, suP_correct);

% Plot errors
figure(17); surf(errors_ap_W); title('ap\_W Error'); colorbar;
figure(18); surf(errors_ap_N); title('ap\_N Error'); colorbar;
figure(19); surf(errors_ap_E); title('ap\_E Error'); colorbar;
figure(20); surf(errors_ap_S); title('ap\_S Error'); colorbar;
figure(21); surf(errors_ap_P); title('ap\_P Error'); colorbar;
figure(22); surf(errors_suP); title('suP Error'); colorbar;

%6.- SOLVE PRESSURE CORRECTION

% Compute errors
errors_p_prime = matrix_err_percentage(p_prime, p_prime_correct);
errors_rsid_p = matrix_err_percentage(rsid_p, rsid_p_correct);
errors_err_p = matrix_err_percentage(err_p, err_p_correct);

% Plot errors
figure(23); surf(errors_p_prime); title('p\_prime Error'); colorbar;
figure(24); surf(errors_rsid_p); title('rsid\_p Error'); colorbar;
figure(25); plot(errors_err_p); title('err\_p Error');

%7.- CORRECT PRESSURE
%compute errors
errors_p_star = matrix_err_percentage(p_star, p_star_correct);
%plot errors
figure(26); surf(errors_p_star); title('Step 7: p\_star Error'); colorbar;

%8.- CORRECT CELL CENTER VELOCITY

% Compute errors
errors_u_star = matrix_err_percentage(u_star, u_star_correct);
errors_v_star = matrix_err_percentage(v_star, v_star_correct);

figure(27); surf(errors_u_star); title('u\_star Error'); colorbar;
figure(28); surf(errors_v_star); title('v\_star Error'); colorbar;

%9.- CORRECT FACE VELOCITY

% Compute errors
errors_u_face_corr = matrix_err_percentage(u_face, u_face_correct);
errors_v_face_corr = matrix_err_percentage(v_face, v_face_correct);

figure(29); surf(errors_u_face_corr); title('u\_face (corrected) Error'); colorbar;
figure(30); surf(errors_v_face_corr); title('v\_face (corrected) Error'); colorbar;

%10.- P=P_STAR

% Compute errors
errors_p = matrix_err_percentage(p, p_correct);
figure(31); surf(errors_p); title('Step 10: Final Pressure Error'); colorbar;
%}
