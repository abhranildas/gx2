function [mu,v]=gx2stat(w,k,lambda,m,s)
	
	parser = inputParser;
	addRequired(parser,'w',@(x) isreal(x) && isrow(x));
	addRequired(parser,'k',@(x) isreal(x) && isrow(x));
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
	addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
	parse(parser,w,k,lambda,m,s);
	
	mu=dot(w,k+lambda)+m;
	v=2*dot(w.^2,k+2*lambda)+s^2;