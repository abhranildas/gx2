function quad_s=standard_quad(quad,mu,v)
% standardize quadratic coefficients
quad_s.q2=sqrtm(v)*quad.q2*sqrtm(v);
quad_s.q1=sqrtm(v)*(2*quad.q2*mu+quad.q1);
quad_s.q0=mu'*quad.q2*mu+quad.q1'*mu+quad.q0;