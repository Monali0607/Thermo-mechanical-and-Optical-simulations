function [I0,u] = bpm_solver(nx,nz,dx,dz,k0,n0,n,E1)

%%%% Monali Suar %%%%%
% Input for the functions are 
% nx - discretized points along x-axis
% nz - discretized points along z-axis 
% dx - discretized step sixe along x-axis
% dz - discretized step size along z-axis
% ko - wave number 
% n0 - reference refractive index value
% n - spatially dependent refractive index(x,z)
% E1 - Input Gaussian beam (1D along x axis)

%% Output
% I0 - Intensity profile calculation by beam propagation method
% u- field value
u = zeros(nx,nz);
u(1:end,1)= E1;

%PARAMETERS 
alpha=0.5;                                                %scheme parameter
a= ones(1,nx).*(-alpha/(dx)^2);
c= ones(1,nx).*(-alpha/(dx)^2);
r_total= zeros(nx,nz);

for m=1:1:nz-1
    for j= 1:1:nx
        b(j)= (2*1i*k0*n0)/dz + (2*alpha)/(dx)^2 - alpha*(n(j,m+1).^2-n0.^2)*k0^2;
             if (j>1 && j<nx) 
                r(j)= (1-alpha)/(dx)^2 *[u(j-1,m)+ u(j+1,m)]+[(1-alpha)*(n(j,m)^2-n0^2)*k0^2-2*(1-alpha)/(dx)^2+2*1i*k0*n0/dz].*u(j,m);
             else  
                r(1) = (1-alpha)/(dx)^2 *[u(2,m)]+[(1-alpha)*(n(1,m)^2-n0^2)*k0^2-2*(1-alpha)/(dx)^2+2*1i*k0*n0/dz].* u(1,m);
    % r(1)= eps;
                r(nx)= (1-alpha)/(dx)^2 *[u(nx-1,m)]+[(1-alpha)*(n(nx,m)^2-n0^2)*k0^2-2*(1-alpha)/(dx)^2+2*1i*k0*n0/dz].* u(nx,m);
    % r(NX) = eps ; 
             end    
    end
    r_total(:,m)= r;
    u(:,m+1)=tridiag(b,a,c,r_total(:,m)); 
end

% Transperant Boundary Condition
 for m = 1:1:nz-1
    boundary_1 = u(nx,m)/u(nx-1,m);        % for upper boundary
    k1 = log(boundary_1)*(-1i)/dx ;
        if real(k1) >= 0
             u(nx, m+1) = u(nx-1,m+1)* exp(1*1i*k1*dx);
        else real(k1)< 0;
             k1 = 0.0;
             u(nx, m+1) = u(nx-1,m+1)*exp(1i*k1*dx);
        end

   boundary_2= u(1,m)/u(2,m);              % for Lower boundary
   k2 = log(boundary_2)*(-1i)/dx ;
     if real(k2) >= 0
         u(1,m+1) = u(2,m+1)* exp(1*1i*k2*dx);
     else real(k2)<=0
         k2 = 0.0;
         u(1, m+1) = u(2,m+1)*exp(1i*k2*dx);
     end
 end
I0= (abs(u)).^2 ; 
end