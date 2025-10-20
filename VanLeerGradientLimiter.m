function [psi] = VanLeerGradientLimiter(r)
%Van Albada Gradient limiter for MUSCL SCHEME
    psi=(r +abs(r))/(1+ abs(r));
end