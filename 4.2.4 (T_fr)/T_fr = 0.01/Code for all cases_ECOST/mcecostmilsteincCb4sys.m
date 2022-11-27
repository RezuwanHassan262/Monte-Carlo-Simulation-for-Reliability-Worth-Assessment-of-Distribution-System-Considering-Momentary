%
% Test weak convergence of Milestein method

%

tic
close all;clear all;
% 
% Monte Carlo simulation 
% 
S0f39  = 1/0.096;   % Failure rate*repair time of 0.6km
S0f4875  = 1/0.12;   % Failure rate*repair time of 0.75km
S0f52  = 1/0.128;   % Failure rate*repair time of 0.8km

S0ft  = 1/0.01;   % Failure rate of x_former

S0r1  =0.5; %switching time

S0r5  =3.0; %repair time of line

S0rt  =10.0; %repair time of x_former


r   = 0.01;  % 
sig = 0.8;   % 


a391=S0r1/S0f39;%component no. 1 & repair 1 hr
a395=S0r5/S0f39; %component no. 1 & repair 5 hr
a48751=S0r1/S0f4875; a48755=S0r5/S0f4875;
a521=S0r1/S0f52; a525=S0r5/S0f52;
at=S0rt/S0ft; 

T   = 1;
%N2  = 19000;   % original total number of Monte Carlo paths % change according to required (228*250/3)
%N2  = 570000; 
N2 = 1000;   % number of paths at a time 
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

ep1=	0.001490572	;
ep2=	0.001489533	;
ep3=	0.001490572	;
ep4=	0.001352957	;
ep5=	0.001344156	;
ep6=	0.057679497	;
ep7=	0.057308583	;
ep8=	0.114702759	;
ep9=	0.114702759	;
ep10=	0.114702759	;
ep11=	0.001490218	;
ep12=	0.001490572	;
ep13=	0.001490572	;
ep14=	0.001486792	;
ep15=	0.001486792	;
ep16=	0.063319514	;
ep17=	0.063319514	;
ep18=	0.001488878	;
ep19=	0.001489872	;
ep20=	0.001488878	;
ep21=	0.001488878	;
ep22=	0.001489872	;
ep23=	0.001488878	;
ep24=	0.063407394	;
ep25=	0.063449277	;
ep26=	0.114702759	;
ep27=	0.114702759	;
ep28=	0.114702759	;
ep29=	0.114702759	;
ep30=	0.114702759	;
ep31=	0.114702759	;
ep32=	0.001489872	;
ep33=	0.001489872	;
ep34=	0.001491305	;
ep35=	0.001486249	;
ep36=	0.001486249	;
ep37=	0.001486249	;
ep38=	0.063296637	;

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
       

 Ff4875 =N2/sum(-Sf4875.*(log(U4875))); %failure of component 0.75
 Ff39 =N2/sum(-Sf39.*(log(U39))); %failure of component 0.60
 Ff52 =N2/sum(-Sf52.*(log(U52))); %failure of component 0.80
  
  Fft =N2/sum(-Sft.*(log(Utf))); %failure of component x former
  
  
  Ffp1=Ff4875+Ff52+Ff52+Ff4875+Ff39+Ff39+Fft;
  Ffp2=Ff4875+Ff52+Ff52+Ff4875+Ff39+Ff4875+Fft;
  Ffp3=Ff4875+Ff52+Ff52+Ff4875+Ff39+Ff39+Fft;
  Ffp4=Ff4875+Ff52+Ff52+Ff4875+Ff39+Ff52+Fft;
  Ffp5=Ff4875+Ff52+Ff52+Ff4875+Ff39+Ff4875+Fft;
  Ffp6=Ff4875+Ff52+Ff52+Ff4875+Ff39+Ff52+Fft;
  Ffp7=Ff4875+Ff52+Ff52+Ff4875+Ff39+Ff4875+Fft;
  
  Ffp8=Ff52+Ff52+Ff39+Ff39;
  Ffp9=Ff52+Ff52+Ff39+Ff4875; 
  Ffp10=Ff52+Ff52+Ff39+Ff52;
  
  Ffp11=Ff4875+Ff39+Ff52+Ff52+Ff39+Ff52+Fft;
  Ffp12=Ff4875+Ff39+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp13=Ff4875+Ff39+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp14=Ff4875+Ff39+Ff52+Ff52+Ff39+Ff39+Fft;
  Ffp15=Ff4875+Ff39+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp16=Ff4875+Ff39+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp17=Ff4875+Ff39+Ff52+Ff52+Ff39+Ff39+Fft;
 
  Ffp18=Ff52+Ff52+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp19=Ff52+Ff52+Ff52+Ff52+Ff39+Ff39+Fft;
  Ffp20=Ff52+Ff52+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp21=Ff52+Ff52+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp22=Ff52+Ff52+Ff52+Ff52+Ff39+Ff39+Fft;
  Ffp23=Ff52+Ff52+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp24=Ff52+Ff52+Ff52+Ff52+Ff39+Ff4875+Fft;
  Ffp25=Ff52+Ff52+Ff52+Ff52+Ff39+Ff39+Fft;
  
  Ffp26=Ff52+Ff39+Ff4875+Ff4875;
  Ffp27=Ff52+Ff39+Ff4875+Ff52;
  Ffp28=Ff52+Ff39+Ff4875+Ff39;
  
  Ffp29=Ff4875+Ff52+Ff52+Ff39;
  Ffp30=Ff4875+Ff52+Ff52+Ff4875;
  Ffp31=Ff4875+Ff52+Ff52+Ff39;
  
  Ffp32=Ff4875+Ff39+Ff4875+Ff4875+Ff4875+Ff52+Fft;
  Ffp33=Ff4875+Ff39+Ff4875+Ff4875+Ff4875+Ff52+Fft;
  Ffp34=Ff4875+Ff39+Ff4875+Ff4875+Ff4875+Ff39+Fft;
  Ffp35=Ff4875+Ff39+Ff4875+Ff4875+Ff4875+Ff52+Fft;
  Ffp36=Ff4875+Ff39+Ff4875+Ff4875+Ff4875+Ff39+Fft;
  Ffp37=Ff4875+Ff39+Ff4875+Ff4875+Ff4875+Ff52+Fft;
  Ffp38=Ff4875+Ff39+Ff4875+Ff4875+Ff4875+Ff39+Fft;
  
  P=((Ffp1*Cp1*ep1)+(Ffp2*Cp2*ep2)+(Ffp3*Cp3*ep3)+(Ffp4*Cp4*ep4)+(Ffp5*Cp5*ep5)+(Ffp6*Cp6*ep6)+(Ffp7*Cp7*ep7)+(Ffp8*Cp8*ep8)+(Ffp9*Cp9*ep9)...
      +(Ffp10*Cp10*ep10)+(Ffp11*Cp11*ep11)+(Ffp12*Cp12*ep12)+(Ffp13*Cp13*ep13)+(Ffp14*Cp14*ep14)+(Ffp15*Cp15*ep15)+(Ffp16*Cp16*ep16)+(Ffp17*Cp17*ep17)...
      +(Ffp18*Cp18*ep18)+(Ffp19*Cp19*ep19)+(Ffp20*Cp20*ep20)+(Ffp21*Cp21*ep21)+(Ffp22*Cp22*ep22)+(Ffp23*Cp23*ep23)+(Ffp24*Cp24*ep24)+(Ffp25*Cp25*ep25)...
      +(Ffp26*Cp26*ep26)+(Ffp27*Cp27*ep27)+(Ffp28*Cp28*ep28)+(Ffp29*Cp29*ep29)+(Ffp30*Cp30*ep30)+(Ffp31*Cp31*ep31)+(Ffp32*Cp32*ep32)...
      +(Ffp33*Cp33*ep33)+(Ffp34*Cp34*ep34)+(Ffp35*Cp35*ep35)+(Ffp36*Cp36*ep36)+(Ffp37*Cp37*ep37)+(Ffp38*Cp38*ep38));
  
  
   %   
      
      
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
% disp (Ffp1); disp (Ffp2);disp (Ffp3);disp (Ffp4);disp (Ffp5);disp (Ffp6);disp (Ffp7);disp (Ffp8);disp (Ffp9);disp (Ffp10);disp (Ffp11);disp (Ffp12);
% disp (Ffp13);disp (Ffp14);disp (Ffp15);disp (Ffp16);disp (Ffp17);disp (Ffp18);disp (Ffp19);disp (Ffp20);disp (Ffp21);disp (Ffp22);
 
 toc