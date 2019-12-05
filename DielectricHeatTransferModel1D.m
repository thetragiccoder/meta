%Initialise Variables:

tNew=60 %seconds, simulation time
dt=1 %seconds, timestep

L=0.05; %metres, slab thickness 
nx=20; %metres, number of steps in x direction
dx=L/nx;%metres, size of step in x direction
x=dx/2:dx:L-dx/2;%metres, position in x direction

k=6 %thermal conductivity
rho=6020 %density
C=527 %specific heat capacity
const=(k)/(rho*C) %thermal diffusivity

T0=298; %Kelvin, initial temperature of slab
TW=500; %Kelvin, initial tempature of left surface
TE=298; %Kelvin, initial tempature of right surface
TC=408.15; %Kelvin, first order transition temperature

K=1.6*(10^5)%Curie constant
g=-1
G=1

E0 = 500

%Create Matrix with Initial values:
T=ones(nx,1)*T0;

TNew=zeros(nx,1);

E=ones(nx,1)*E0

ENew=zeros(nx,1)
   ENew(1) = 500
    ENew(nx) = 500
    
t=0:dt:tNew

for i = 1:length(t)
    for j = 2:nx-1
        TNew(j) = const*((T(j+1)-T(j))/dx^2 + (T(j-1)-T(j))/dx^2);% Master equation
         
         if (TNew>TC)
    ENew(j) = K/(T(j)-TC)
         else 
         ENew(j) = K/(2*(TC-T(j)))
          %ENew(j) = (4*G/(3*g.^2)) + (K/(8*TNew(j)-TC))
      end
 TNew(1) = const*((T(2)-T(1))/dx^2 + (TW-T(1))/dx^2); %right boundary condition
 TNew(nx) = const*((TE-T(nx))/dx^2 + (T(nx-1)-T(nx))/dx^2); %left boundary condition
    end
    
 E = E+ENew
 T = T+(TNew*dt);

    
 
end
 

 figure (1) % new figure
[hAx,hLine1,hLine2] = plotyy(x,T,x,E);
    title('Dielectric Constant and Temperature of a Heated 1D Slab of BaTiO3')
    xlabel('Distance (m)')
    ylabel(hAx(1),'Temperature(K)') % left y-axis 
ylabel(hAx(2),'Dielectric Constant(\epsilon)') % right y-axis
    %plot (x,E,'r','Linewidth',2)
    %pause(1)
    
    figure(2)
    plot(T,E)
    title('Dielectric Constant as a Function of Temperature of BaTiO3')
    xlabel('Temperature (K)')
    ylabel('Dielectric Constant (\epsilon)')

