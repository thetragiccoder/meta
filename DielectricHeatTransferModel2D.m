%Initialise Variables:

tNew=60 %seconds, simulation time
dt=1 %seconds, timestep

Lx=500; %metres, slab thickness in x direction
Ly=500; %metres, slab thickness in y direction
nx=20; %metres, number of steps in x direction
ny=20; %metres, number of steps in y direction
dx=Lx/nx;%metres, size of step in x direction
dy=Ly/ny;%metres, size of step in y direction
x=dx/2:dx:Lx-dx/2;%metres, position in x direction
y=dy/2:dy:Ly-dy/2;%metres, position in y direction
[X,Y]=meshgrid(x,y);%create mesh for 3D plot

k=6 %thermal conductivity
rho=6020 %density
C=527 %specific heat capacity
const=(k)/(rho*C) %thermal diffusivity

T0=298; %Kelvin, initial temperature of slab
TW=500; %Kelvin, initial tempature of left surface
TE=298; %Kelvin, initial tempature of right surface
TN=500; %Kelvin, initial tempature of top surface
TS=298; %Kelvin, initial tempature of bottom surface
TC=ones(nx,ny)*408.15; %Kelvin, first order transition temperature

K=1.6*(10^5)%Curie constant
g=-1
G=1

E0 = 500

%Create Matrix with Initial values:
T=ones(nx,ny)*T0;

TNew=zeros(nx,ny);

E=ones(nx,ny)*E0

ENew=zeros(nx,ny)
   ENew(nx) = 500
    ENew(ny) = 500
    
t=0:dt:tNew

for n = 1:length(t)
    for i = 2:nx-1
    for j = 2:nx-1
        TNew(i,j) = const*((T(j+1)-T(j))/dx^2 + (T(j-1)-T(j))/dx^2 +...
            (T(i+1)-T(i))/dx^2 + ((T(i-1)-T(i))/dx^2));% Master equation
         
       if (TNew(i,j)>TC(i,j))
    ENew(i,j) = K./((T(i,j)-TC(i,j)))
        else 
         ENew(i,j) = K./(2.*(TC(i,j)-T(i,j)))
          %ENew(i,j) = (4*G/(3*g.^2)) + (K./(8.*TNew(i,j)-TC(i,j)))
        end
           
    end
        
    end
    
     TNew(1,:) = const*((T(2,:)-T(1,:))/dx^2 + (TW-T(1,:))/dx^2); %right boundary condition
    TNew(nx,:) = const*((TE-T(nx,:))/dx^2 + (T(nx-1,:)-T(nx,:))/dx^2); %left boundary condition
    TNew(:,1) = const*((T(:,2)-T(:,1))/dx^2 + (TN-T(:,1))/dx^2); %top boundary condition
    TNew(:,ny) = const*((TS-T(nx,:))/dx^2 + (T(nx-1,:)-T(nx,:))/dx^2); %bottom boundary condition
     E = E+ENew
    T = T+(TNew*dt);
    
    figure(1) 
    mesh (x,y,T)
     pause(1)
     
     figure(2)
     mesh (x,y,E)
     pause(1)

    end
    
   
  
   %plot result
    %figure(1)
    %plot(T,E)
   % plot(x,T,'b','Linewidth', 2)
    %title('Dielectric Constant of Heated 1D Slab of BaTiO3')
    %xlabel('Distance (m)')
    %ylabel('Temperature (K)')
    %hold on
    %yyaxis right
    %ylabel('Dielectric constant')
    %plot (x,E,'r','Linewidth',2)
    %pause(1)
 %plot(x,T,'b','Linewidth', 2)
    %title('Dielectric Constant of Heated 1D Slab of BaTiO3')
    %xlabel('Distance (m)')
    %ylabel('Temperature (K)')
    %hold on
    %yyaxis right
    %ylabel('Dielectric constant')
    %plot (x,E,'r','Linewidth',2)
    
    
