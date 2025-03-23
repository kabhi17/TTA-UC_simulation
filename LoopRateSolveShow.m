% MATALB script for transient kinetic analysis developed by Abhishek
% Kalpattu at the Fourkas Lab, UMD
%% %% Variable definition
% Define all variables here
kfl= % Annihilator fluorescence rate constant
kNR= % Non=radiative annihilator singlet quenching rate constant
ktt= % TTA rate constant
ksens= % TET rate constant
krsens= % Reverse TET rate constant
kttS= % Sensitizer TTA rate constant
kts= % Senstiizer intrinsic triplet quenching rate constant
kt= % Annihilator intrinsic tript quenching rate constant
kqS= % Sensitizer self-quenching rate constant
kq= % Dynamic quenching rate constant
kfret= % FRET rate constant
kic= % Internal conversion rate constant
kph= % Sensitizer phosphorescence rate constant
kext= % Parasitic sensitizer rexcitation rate constant
krisc= % RISC rate constant in the annihilator
kex= % Sensitizer excitation rate constant 
S0= % Sensitizer conc.
A0= % Annihilator conc.
dt= % Timestep of the simulation
%% Time-resolved PL simulation 
Sto=[]; % Array for storing sensitizer phosphorescence values w.r.t time
Sto1=[]; % Array for storing upconverted fluorescence values w.r.t time
% These are arrays for storing information on the conc. of various excited
% state species
Fro=[];
Kro=[];
Gro=[];
% Initiatlizing the conc. of excited state species
As=0;
Ss=0;
Aa=0;
Asa=0;
AaC=[];
AsC=[];
SsC=[];
AsaC=[];
I=2e3; % This is the irradiance of the excitation source as it is 'turned on'
for i=0:100000 % For-loop block 1
SsC=[SsC,Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt];
Ss=Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt;
AsC=[AsC,As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Asa-ksens*Ss*Aa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt];
As=As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Aa-ksens*Ss*Asa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt;
AsaC=[AsaC,Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt];
Asa=Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt;
AaC=[AaC,Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt];
Aa=Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt;
end
Sto=[Sto,kph*SsC];Sto1=[Sto1,kfl*AaC];Fro=[Fro,SsC];Kro=[Kro,AsC];Gro=[Gro,AsaC];
% Arrays reinitialized but the last value for [3A*], [1A*], and [3S*] are
% kept
AaC=[];
AsC=[];
SsC=[];
AsaC=[];
for i=0:100000 % For-loop block 2 (Irradiance is unchanged at preset value)
SsC=[SsC,Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt];
Ss=Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt;
AsC=[AsC,As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Asa-ksens*Ss*Aa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt];
As=As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Aa-ksens*Ss*Asa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt;
AsaC=[AsaC,Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt];
Asa=Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt;
AaC=[AaC,Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt];
Aa=Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt;
end
Sto=[Sto,kph*SsC];Sto1=[Sto1,kfl*AaC];Fro=[Fro,SsC];Kro=[Kro,AsC];Gro=[Gro,AsaC];
I=0; % Irradiance is set to 0, simulating the excitation beam being turned off
AaC=[];
AsC=[];
SsC=[];
AsaC=[];
for i=0:100000 % For-loop block 3
SsC=[SsC,Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt];
Ss=Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt;
AsC=[AsC,As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Asa-ksens*Ss*Aa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt];
As=As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Aa-ksens*Ss*Asa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt;
AsaC=[AsaC,Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt];
Asa=Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt;
AaC=[AaC,Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt];
Aa=Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt;
end
Sto=[Sto,kph*SsC];Sto1=[Sto1,kfl*AaC];Fro=[Fro,SsC];Kro=[Kro,AsC];Gro=[Gro,AsaC];
AaC=[];
AsC=[];
SsC=[];
AsaC=[];
for i=0:100000 % For-loop block 4
SsC=[SsC,Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt];
Ss=Ss+(kex*I*S0+krsens*S0*As+kfret*S0*Aa-(Ss*(kex*I+krsens*As+kfret*Aa+ksens*(A0-As-Aa-Asa)+kts+kph+kqS*(S0-Ss)+kttS*Ss)))*dt;
AsC=[AsC,As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Asa-ksens*Ss*Aa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt];
As=As+(ksens*Ss*A0-krsens*As*(S0-Ss)-ksens*Ss*As-ksens*Ss*Aa-ksens*Ss*Asa-2*ktt*(As^2)+kic*Asa-kt*As-kq*(S0-Ss)*As)*dt;
AsaC=[AsaC,Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt];
Asa=Asa+(0.75*ktt*(As^2)-kic*Asa-krisc*Asa)*dt;
AaC=[AaC,Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt];
Aa=Aa+(0.25*ktt*(As^2)+krisc*Asa+kfret*Ss*Aa-kfret*S0*Aa-(kfl+kNR)*Aa)*dt;
end
Sto=[Sto,kph*SsC];Sto1=[Sto1,kfl*AaC];Fro=[Fro,SsC];Kro=[Kro,AsC];Gro=[Gro,AsaC];
t=(0:dt:length(Sto)*dt-dt); % Total time length of the simulation
%% Data plotting
hold on; plot(t-2e-4,log10(Sto./max(Sto)),'LineWidth',1.5); % Plot normalized sensitizer phosphorescence as a function of time
hold on; plot(t-2e-4,log10(Sto1./max(Sto1)),'LineWidth',1.5); % Plot normalized upconverted fluorescence as a function of time
% Calculating of triplet lifetimes 
LifeT=[LifeT,1./vpa((log10(Sto1(400003)./max(Sto1))-log10(Sto1(400004)./max(Sto1)))./(dt))]; % Annihilator triplet lifetime from the UCPL decay tail
LifeTP=[LifeTP,1./vpa((log10(Sto(400003)./max(Sto))-log10(Sto(400004)./max(Sto)))./(dt))]; % Effective annihilator triplet lifetime from the sensitizer phosphorescence decay tail

