function F=nonlinearsolver(theta,PI,Z)
global A kappa beta alpha mu b delta; 
global PI Z;
F=(kappa/(beta*A)).*theta.^(alpha)-PI*((1-mu).*(Z-b)-(kappa*mu).*theta+((1-delta)*kappa/A).*theta.^(alpha));
end