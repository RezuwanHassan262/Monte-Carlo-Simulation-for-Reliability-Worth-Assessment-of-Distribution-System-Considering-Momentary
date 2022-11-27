%
% Test weak convergence of Milestein method

%
close all;

clear all;

tic

% 
% Monte Carlo simulation 
% 

S0f39  = 1/0.045;   % Failure rate*length of 0.6km
S0f4875  = 1/0.05625;   % Failure rate*length of 0.75km
S0f52  = 1/0.06;   % Failure rate*length of 0.8km

S0ft  = 1/0.1;   % Failure rate of x_former

S0r1  =0.5; %switching time

S0r5  =3.0; %repair time of line

S0rt  =10.0; %repair time of x_former

r   = 0.02;  % 
sig = 0.1;   %

 

a391=S0r1/S0f39;%component no. 1 & repair 1 hr
a395=S0r5/S0f39; %component no. 1 & repair 5 hr
a48751=S0r1/S0f4875; a48755=S0r5/S0f4875;
a521=S0r1/S0f52; a525=S0r5/S0f52;
at=S0rt/S0ft;

T   = 1;
N2  =184000;   % total number of Monte Carlo paths % change according to required 
%M2 = 1000;   % number of paths at a time 
Eps   = 0.03;% original %Change according to accuracy required
%Eps   = 0.01;
%Eps   = 0.05;
%Eps   = 0.001;
L =3;% change according to MLMC
Ns = 2^L;
hf = T/Ns;

Cp1=0.545; Cp2=0.545;Cp3=0.545;Cp4=0.545; Cp5=0.5;Cp6=0.415; Cp7=0.415;Cp8=1;Cp9=1.5;Cp10=1;
Cp11=0.545; Cp12=0.545;Cp13=0.545;Cp14=0.5; Cp15=0.5;Cp16=0.415; Cp17=0.415;Cp18=0.545;Cp19=0.545;Cp20=0.545;
Cp21=0.545;Cp22=0.5;Cp23=0.5;Cp24=0.415; Cp25=0.415;Cp26=1; Cp27=1;Cp28=1;Cp29=1;Cp30=1;
Cp31=1.5;Cp32=0.545;Cp33=0.545;Cp34=0.545; Cp35=0.545;Cp36=0.5; Cp37=0.5;Cp38=0.415;


  sum1 = 0;
  sum2 = 0;

  for i = 1:N2
  %  m2 = min(M2,M-m+1); 

 Sf39 = S0f39*ones(1,N2); 
  Sf4875 = S0f4875*ones(1,N2); 
  Sf52 = S0f52*ones(1,N2);
    
  Sft = S0ft*ones(1,N2);
   
   
 U4875  = rand(1,N2); U52  = rand(1,N2);U39  = rand(1,N2);
  Utf  = rand(1,N2);Utr  = rand(1,N2);U1r  = rand(1,N2);U5r  = rand(1,N2);
  
if U4875==1;
   U4875=0.999;
end
if U52==1;
   U52=0.999;
end
if U39==1;
   U39=0.999;
end

if Utf==1;
   Utf=0.999;
end
if Utr==1;
   Utr=0.999;
end
if U1r==1;
   U1r=0.999;
end
if U5r==1;
   U5r=0.999;
end
    for n = 1:Ns
      dWf = sqrt(hf)*randn(1,N2);
 Sf39  = Sf39 + r*Sf39*hf + sig*Sf39.*dWf ...
        + 0.5*sig^2*Sf39.*(dWf.^2-hf);%Component 0.6
    
   Sf4875  = Sf4875 + r*Sf4875*hf + sig*Sf4875.*dWf ...
       + 0.5*sig^2*Sf4875.*(dWf.^2-hf);%Component 0.75 
  
   Sf52  = Sf52 + r*Sf52*hf + sig*Sf52.*dWf ...
       + 0.5*sig^2*Sf52.*(dWf.^2-hf);%Component 0.8 
   
     Sft  = Sft + r*Sft*hf + sig*Sft.*dWf ...
       + 0.5*sig^2*Sft.*(dWf.^2-hf);%Component x former
    end
       
 Uf391=sum(-a391.*Sf39.*(log(U1r)))/sum((-Sf39.*(log(U39)))+(-a391.*Sf39.*(log(U1r)))/8760);%unavailability of component 1 repair 1 hour
  Uf395=sum(-a395.*Sf39.*(log(U5r)))/sum((-Sf39.*(log(U39)))+(-a395.*Sf39.*(log(U5r)))/8760);%unavailability of component 1 repair 5 hour
  Uf48751=sum(-a48751.*Sf4875.*(log(U1r)))/sum((-Sf4875.*(log(U4875)))+(-a48751.*Sf4875.*(log(U1r)))/8760);%unavailability of component 1 repair 1 hour
  Uf48755=sum(-a48755.*Sf4875.*(log(U5r)))/sum((-Sf4875.*(log(U4875)))+(-a48755.*Sf4875.*(log(U5r)))/8760);%unavailability of component 1 repair 5 hour
  Uf521=sum(-a521.*Sf52.*(log(U1r)))/sum((-Sf52.*(log(U52)))+(-a521.*Sf52.*(log(U1r)))/8760);%unavailability of component 1 repair 1 hour
  Uf525=sum(-a525.*Sf52.*(log(U5r)))/sum((-Sf52.*(log(U52)))+(-a525.*Sf52.*(log(U5r)))/8760);%unavailability of component 1 repair 5 hour
  
  Uft=sum(-at.*Sft.*(log(Utr)))/sum((-Sft.*(log(Utf)))+(-at.*Sft.*(log(Utr)))/8760);%unavailability of component 1 repair 5 hour

  Ufp1=Uf48751+Uf521+Uf521+Uf48751+Uf391+Uf395+Uft;
  Ufp2=Uf48751+Uf521+Uf521+Uf48751+Uf391+Uf48755+Uft;
  Ufp3=Uf48751+Uf521+Uf521+Uf48751+Uf391+Uf395+Uft;
  Ufp4=Uf48751+Uf521+Uf521+Uf48751+Uf391+Uf525+Uft;
  Ufp5=Uf48751+Uf521+Uf521+Uf48751+Uf391+Uf48755+Uft;
  Ufp6=Uf48751+Uf521+Uf521+Uf48751+Uf391+Uf525+Uft;
  Ufp7=Uf48751+Uf521+Uf521+Uf48751+Uf391+Uf48755+Uft;
  
  Ufp8=Uf521+Uf521+Uf391+Uf395;
  Ufp9=Uf521+Uf521+Uf391+Uf48755; 
  Ufp10=Uf521+Uf521+Uf391+Uf525;
  
  Ufp11=Uf48751+Uf391+Uf521+Uf521+Uf391+Uf525+Uft;
  Ufp12=Uf48751+Uf391+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp13=Uf48751+Uf391+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp14=Uf48751+Uf391+Uf521+Uf521+Uf391+Uf395+Uft;
  Ufp15=Uf48751+Uf391+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp16=Uf48751+Uf391+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp17=Uf48751+Uf391+Uf521+Uf521+Uf391+Uf395+Uft;
  
  Ufp18=Uf521+Uf521+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp19=Uf521+Uf521+Uf521+Uf521+Uf391+Uf395+Uft;
  Ufp20=Uf521+Uf521+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp21=Uf521+Uf521+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp22=Uf521+Uf521+Uf521+Uf521+Uf391+Uf395+Uft;
  Ufp23=Uf521+Uf521+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp24=Uf521+Uf521+Uf521+Uf521+Uf391+Uf48755+Uft;
  Ufp25=Uf521+Uf521+Uf521+Uf521+Uf391+Uf395+Uft;
  
  Ufp26=Uf521+Uf391+Uf48751+Uf48755;
  Ufp27=Uf521+Uf391+Uf48751+Uf525;
  Ufp28=Uf521+Uf391+Uf48751+Uf395;
  
  Ufp29=Uf48751+Uf521+Uf521+Uf395;
  Ufp30=Uf48751+Uf521+Uf521+Uf48755;
  Ufp31=Uf48751+Uf521+Uf521+Uf395;
  
  Ufp32=Uf48751+Uf391+Uf48751+Uf48751+Uf48751+Uf525+Uft;
  Ufp33=Uf48751+Uf391+Uf48751+Uf48751+Uf48751+Uf525+Uft;
  Ufp34=Uf48751+Uf391+Uf48751+Uf48751+Uf48751+Uf395+Uft;
  Ufp35=Uf48751+Uf391+Uf48751+Uf48751+Uf48751+Uf525+Uft;
  Ufp36=Uf48751+Uf391+Uf48751+Uf48751+Uf48751+Uf395+Uft;
  Ufp37=Uf48751+Uf391+Uf48751+Uf48751+Uf48751+Uf525+Uft;
  Ufp38=Uf48751+Uf391+Uf48751+Uf48751+Uf48751+Uf395+Uft;
  
   P=((Ufp1*Cp1)+(Ufp2*Cp2)+(Ufp3*Cp3)+(Ufp4*Cp4)+(Ufp5*Cp5)+(Ufp11*Cp11)+(Ufp12*Cp12)+(Ufp13*Cp13)+(Ufp14*Cp14)...
      +(Ufp15*Cp15)+(Ufp18*Cp18)+(Ufp19*Cp19)+(Ufp20*Cp20)+(Ufp21*Cp21)+(Ufp22*Cp22)...
      +(Ufp23*Cp23)+(Ufp32*Cp32)+(Ufp33*Cp33)+(Ufp34*Cp34)+(Ufp35*Cp35)+(Ufp36*Cp36)+(Ufp37*Cp37));
  
%disp(P);
sum1 = sum1 + sum(P);
%disp(sum1/i);
sum2 = sum2 + sum(P.^2); 
a=sum(i);
  ml  = sum1/a; %MEAN
  Vl  = sum2/a - ml^2;%VARIANCE
  e1=Vl/a;
  if e1>=4.400e-004 && e1<=4.50000e-004   %0.03
  %if e1>=4.900e-005 && e1<=5.00000e-005   %0.01
  %if e1>=1.200e-003 && e1<=1.25000e-003   %0.05
  %if e1>=4.900e-007 && e1<=5.0000e-007   %0.001   
  break 
  end
  end
 % ml  = sum1/N2; %MEAN
 % Vl  = sum2/N2 - ml^2;%VARIANCE
 e=(Eps^2)/2;
 % e1=Vl/N2;
  disp(ml);
  disp(Vl);
  disp(P);
  
  disp(e);
  disp(e1);
  if e1<=e
 disp('Ok');
 disp(a);
  end 
 %end
 toc