function [psi] = VanAlbadaGradientLimiter(r)
%Van Albada Gradient limiter for MUSCL SCHEME
    psi=(r +r^2)/(1+ r^2);
end