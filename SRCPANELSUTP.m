% This code is intended to simulate 2-D potential flows over nonlifting 
% bodies.  Thefoundations of this code can be found on John Anderson's book
% Fundamentals of Aerodynamics, Fourth Edittion
%
% It basically takes the geometry and covers it with source panels and
% compute the flowfield around it
%

% Input:  The arbitrary shape
% Output: Velocity potential
%         Velocity Vectors
%         Streamlines
%         Cp distribution
% Oscar Garibaldi Ph.D., 2015

clear all
close all
clc
%% LOADING GEOMETRY
filename = 'DIAMANTE.txt'; % Name of the text file containing the body in x and y coordinates. Modify at will!

BODY = load(filename);

%% GENERATION OF PANELS

% Matrix Panel format
%  [ Xj Yj Xj+1 Yj+1 ORIENTATION NORMAL Panel Lenght ]

[n,m] = size(BODY);

PANEL = zeros(n-1,10);

for i=1:n-1
    
    PHI =180/pi*atan2(BODY(i+1,2) - BODY(i,2),BODY(i+1,1) - BODY(i,1));
    S = sqrt( (BODY(i+1,2) - BODY(i,2))^2 + (BODY(i+1,1) - BODY(i,1))^2 );
    
    xi =( BODY(i,1) + BODY(i+1,1) )/2;
    yi =( BODY(i,2) + BODY(i+1,2))/2;
    
    PANEL(i,:)= [ BODY(i,1) BODY(i,2) BODY(i+1,1) BODY(i+1,2) PHI  PHI + 90 S xi yi 0];   
          
end

 plot(PANEL(:,8),PANEL(:,9),'o');  % Plotting the Control points


%% Computing the integral (Influence of each panel of each control points)

Vinf = 1;

PANELJ = PANEL; 

n = size(PANELJ,1);

IIJ =ones(n,n)*pi;

FIJ =zeros(n,n)*pi;

VCOS=zeros(n,1);
VSIN=zeros(n,1);


 for i = 1:n

     PANELI = [PANELJ(i,:)];
     
    xi =(PANELI(1,1) + PANELI(1,3))/2;
    yi =(PANELI(1,2) + PANELI(1,4))/2;

    PHIi = (PANELI(1,5))*pi/180;

    BETAi = (PANELI(1,6))*pi/180;


    for j=1:n
        if j~=i
            
           
            XJ = PANELJ(j,1);
            XJ1 = PANELJ(j,3);
            YJ = PANELJ(j,2);
            YJ1 = PANELJ(j,4);
            
            PHIJ = PANELJ(j,5)*pi/180;
            BETAJ = PANELJ(j,6)*pi/180;
            SJ = PANELJ(j,7);
            
            
            A = -(xi - XJ)*cos(PHIJ) - (yi - YJ)*sin(PHIJ);
            B = (xi - XJ)^2 + (yi - YJ)^2;
            C = sin(PHIi - PHIJ);
            D = (yi - YJ)*cos(PHIi) - (xi - XJ)*sin(PHIi);
            E = sqrt( B - A^2);
         IIJ(i,j) = (C/2)*log((SJ^2 + 2*A*SJ + B)/B) + (D-A*C)/E*(atan((SJ+A)/E) - atan(A/E)); % This is the integral for lambdas
         
         FIJ(i,j) = (D - A*C)/(2*E)*log((SJ^2 + 2*A*SJ + B)/B) - C*(atan((SJ+A)/E) - atan(A/E)); %This is the integral for the velocity
         
        end
    end

    VCOS(i,1) = Vinf*cos(BETAi)*2*pi;
    
    VSIN(i,1) = Vinf*sin(BETAi);

 end
 
 % Encontrando las soluciones de Source Strenghts
 lambda = IIJ\-VCOS;
 
 sum(lambda) % Verification of results, Should be ZERO!

 %% Finding the Tangential Velocities and Cp

 VI=zeros(n,1);

 CPI=zeros(n,1);
 
 for i = 1:n

     VI(i,1) = VSIN(i,1) + sum(lambda'.*FIJ(i,:))/(2*pi);
     
     CPI(i,1) = 1 - (VI(i,1)/Vinf)^2;
     
 end
 
  
