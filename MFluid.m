clear; close all; clc;
%% Input parametrs 
lx = 100; %System length
dx = 0.2; %Spatial grid size
dt = 0.1; %Temporal grid size
ne0 = 1.0; %Equilibrium electron density
ni0 = 1.0; %Equilibrium ion density
nx = floor(lx/dx); %Number of grid points
x0 = 0.5*nx*dx; %Position of the perturbation
wd = 2; %Width of the perturbation
a0 = 0.25; %Amplitude of perturbation (same in this case)
%k_p = 2;
w = 0.70;
tol = 1*exp(-7);
ntime = 1000;
nx4 = nx + 4; %Extra grid points to update the periodic boundary

%% Initializing variables 
x = linspace(-1, double(nx4) - 2, double(nx4))*dx;
x(1)= 0.0;
x(2)=0.0;
x(length(x))=0.0;
x(length(x)-1)=0.0;
%--------------------------------------------------------------------------
%Initializing variables
ne = zeros(nx4,2); %Electron density
ni = zeros(nx4,2); %Ion density
vi = zeros(nx4,2); %Ion velocity
ph = zeros(nx4,2); %Potential
Ex = zeros(nx4,1); %Electric field
PH = zeros(nx4,100); %Electric field
%---------------------------------------------------------------------------
%Perturbation system
noise = a0*rand(nx4,1); %normal random number

%ne(:,1) = ne0 + a0*exp(-1*((x-x0)/wd).^2);
%ne(:,2) = ne0 + a0*exp(-1*((x-x0)/wd).^2);
ne(:,1) = gaussian(x, 2, 10, 50);
ne(:,2) = gaussian(x, 2, 10, 50);

%ni(:,1) = ni0 + a0*exp(-1*((x-x0)/wd).^2);
%ni(:,2) = ni0 + a0*exp(-1*((x-x0)/wd).^2);
ni(:,1) = gaussian(x, 2, 10, 50);
ni(:,2) = gaussian(x, 2, 10, 50);


plot(x(3:nx4 - 2), ni(3:nx4 - 2,1))
%xlim([400,600])
%ylim([-2,2])
%title('Iteration = ' )
%show()
%--------------------------------------------------------------------------

%Computing electrostatic potential when system get perturb
for i = 4 : nx4-2
    ph(i+1,1) = 2*ph(i,1) - ph(i-1,1) - dx*dx*(ni(i,1) - ne(i,1));
end
ph(:,2) = ph(:,1);

jj=1;
for k = 1 : ntime
    %for time 2
    %the continuity equation
    [vi,ni] = continuity_equation(vi,ni,ph,2,dx,dt);
    
    %the poisson equation
    ph(:,2) = poissons_solution(ph(:,2), ni(:,2),dx,w,tol);
    %initial guess
    ph(:,1) = ph(:,2);
    
    %for time 1
    %the continuity equation
    [vi,ni] = continuity_equation(vi,ni,ph,1,dx,dt);
    
    %the poisson equation
    ph(:,1) = poissons_solution(ph(:,1), ni(:,1),dx,w,tol);
    %initial guess
    ph(:,2) = ph(:,1);    
    
 %Plotting or saving the results
    PH(:,jj) = ph(:,1);
    jj = jj + 1;
    if mod(k, 100)==0
            plot(x(3:nx4 - 2), ph(3:nx4 - 2,1))
           % xlim([400,600])
            %ylim([-2,2])
            title(strcat('Iteration = ',int2str(k)))
            %show()
            saveas(gcf, strcat('mfluid_iteration',int2str(k),'.pdf'))
            jj=1;
    end
end   
%-------------------------------------------------------------------------
%Continuity and Momentum Equation
function [x1,y1] = continuity_equation(x,y,z,ic,dx,dt)
%x: Ion velocity
%y: Ion density
%z: Electrostatic potential
%ic: Time index 2 or 1
%dx: Grid size -space
%dt: Grid size -time

%Previous time computation
    jc = 2 - ic;
    fd0 = fd4(x(:,jc+1),dx);
    fd1 = fd4(y(:,jc+1),dx);
    fd2 = fd4(z(:,jc+1),dx);

    for i = 3 : max(size(x))-2
        y(i,ic) = y(i,ic) - dt*(y(i,jc+1)*fd0(i) - x(i,jc+1)*fd1(i));
        x(i,ic) = x(i,ic) - dt*(x(i,jc+1)*fd0(i) + fd2(i));
    end
        
    %periodic boundary 
    x(:,ic)= periodic_boundary(x(:,ic));
    y(:,ic)= periodic_boundary(y(:,ic));
        
    %Filtering numerical noise    
    x(:,ic) = filter3p(x(:,ic));
    y(:,ic) = filter3p(y(:,ic));
    x1 = x;
    y1 = y;
end
   
function [y] = fd4(z, xd)
    y = zeros(1,length(z));
    for i = 3 : (max(size(z))-2)
        y(i) = (8*z(i+1) - 8*z(i-1)-z(i+2)+z(i-2))./(12*xd);
    end
end


function [z1] = filter3p(z)
    for i = 3 : max(size(z))-2
        z(i) = (-1 * z(i - 2) + 4 * z(i - 1) + 10 * z(i)...
                + 4 * z(i + 1) - z(i + 2))./16;
    end
    z = periodic_boundary(z);
    z1= z;
end
%-----------------------------------------------------------------------

%Poisson Equation
function [x1] = poissons_solution(x,z,dx,w,tol)
%x: Electrostatic potential
%z: Ion density
%dx: Grid size-space
%w: SOR method coefficent
%tol: Tolerance for accuracy

%tol = 1*exp(-7);
    pp = 2.0;
    while pp > tol
        xold = x;
        for i = 3 : (max(size(x))-2)
            p0 = 0.5 * (x(i + 1) + x(i - 1) + dx*dx*(z(i) - exp(x(i))));
            x(i) = p0 + w * (p0 - xold(i));
        end
        
    %periodic boundary 
        x = periodic_boundary(x);
        
        diff = xold - x;
        diff = abs(diff);
        pp = max(diff);
        
        x = filter3p(x);
        x1 = x;
    end
end
%---------------------------------------------------------------------------

%Gaussian function
function [y2]= gaussian(x, a, L0, c)
     %It return the Gaussian function output for given array.
     %Parameters
     %------------------------------------------------------
     %x = linspace(1,1000, 1000);
     %n = gaussian(x, 2, 10, 50);
    
     y = zeros(1, length(x));
     for i = 1: max(size(x))
         p = ((x(i) - c) /L0).^2;
         y(i) = a*exp(-1*p);
     end
     y = periodic_boundary(y);
     y2 = y;
 end

%Boundary condition
function [z2] = periodic_boundary(z)
     z(1)= z(length(z)-3);
     z(2)= z(length(z)-2);
     z(length(z)-1)= z(3);
     z(length(z))= z(4);
     z2 = z;
end
