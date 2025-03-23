% MATALB script for simulating the performance of TTA-UC systems, by
% Abhishek Kalpattu at the Fourkas Lab, UMD
%% Variable definition
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
%% General model with all UEL mechanisms
FLA=[]; % Empty initialized matrix for storing steady-state upconverted fluorescence rates
Eff=[]; % Empty initialized matrix for storing upconverted fluorescence quantum yields
EmR=[]; % Empty initialized matrix for storing emission ratios
PH=[]; % Empty initialized matrix for storing steady-state senstiizer phosphorescence rates
for i=-3:0.05:6
    I=10.^i; % Array of irradiances to simulate
    syms Astar Sstar % Symbolic unknown variables (Sstar = sensitizer triplet state conc.), (Astar = annihilator triplet state conc.)
    Aastar=(0.25*ktt*(Astar^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))/(kfret*(S0-Sstar)+kfl+kNR); % Expression for steady-state annihilator single conc.

    Eqn=kex*I*S0-kex*I*Sstar-ksens*Sstar*(A0-Astar-Aastar-0.75*(ktt/(kic+krisc))*(Astar^2))-kqS*(S0-Sstar)*Sstar-kttS*(Sstar.^2)+krsens*(S0-Sstar)*Astar+kfret*Aastar*(S0-Sstar)-(kph+kts)*Sstar==0;
    Sstar=solve(Eqn,Sstar); % Solves for Sstar in terms of Astar

    Aastar=(0.25*ktt*(Astar^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))/(kfret*(S0-Sstar)+kfl+kNR); 
    Eqn2=ksens*Sstar*(A0-Astar-Aastar(1)-0.75*(ktt/(kic+krisc))*(Astar^2))-2*ktt*(Astar^2)+(kic/(krisc+kic))*0.75*ktt*(Astar^2)-kt*Astar-krsens*(S0-Sstar)*Astar-kq*(S0-Sstar)*Astar==0;
    Astar1=vpasolve(Eqn2(1),Astar); % Solves for Astar (Note that there are two equations for Astar that can be solved, only one yields the right answer)

    if Astar1(1)<0 % Logical check to disregard the negative solution for Astar
        Astar2=Astar1(3);
    else 
        Astar2=Astar1(1);
    end

    Sstar=subs(Sstar,Astar,Astar2);
    Astar=Astar2;
    FL=kfl*((0.25*ktt*(Astar.^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))./(kfret*(S0-Sstar)+kfl+kNR)); % Calculated value for upconverted fluorescenc rate
    Aastar=((0.25*ktt*(Astar.^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))./(kfret*(S0-Sstar)+kfl+kNR));
    % The matrices below store all relavant emission data that is to be
    % plotted 
    FLA=[FLA,FL(1)]; % UCPL rate
    PH=[PH,kph*Sstar(1)]; % Sensitizer Phosphorescence
    EmR=[EmR,(kph*Sstar(1))./FL(1)]; % Emission ratios
    Eff=[Eff,FL(1)/(kex*I*(S0-Sstar(1)))]; % QY
    EffEx=[EffEx,FL(1)/I]; % This is the extrinsic qunatum yield
end
%% Parasitic sensitizer reabsorption
FLA=[];
Eff=[];
EmR=[];
PH=[];
EffEx=[];
for i=-3:0.05:6
    I=10.^i;
    syms Astar Sstar
    Aastar=(0.25*ktt*(Astar^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))/(kfret*(S0-Sstar)+kfl+kNR);
    Eqn=kex*I*S0-kex*I*Sstar-ksens*Sstar*(A0-Astar-Aastar-0.75*(ktt/(kic+krisc))*(Astar^2))-kqS*(S0-Sstar)*Sstar-kttS*(Sstar.^2)+krsens*(S0-Sstar)*Astar+kfret*Aastar*(S0-Sstar)-(kph+kts)*Sstar+kexP*(S0-Sstar)*((kfl*Aastar)/(kext+kexP*(S0-Sstar)))==0;
    Sstar=solve(Eqn,Sstar);
    Aastar=(0.25*ktt*(Astar^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))/(kfret*(S0-Sstar)+kfl+kNR);
    Eqn2=ksens*Sstar*(A0-Astar-Aastar(1)-0.75*(ktt/(kic+krisc))*(Astar^2))-2*ktt*(Astar^2)+(kic/(krisc+kic))*0.75*ktt*(Astar^2)-kt*Astar-krsens*(S0-Sstar)*Astar-kq*(S0-Sstar)*Astar==0;
    % Numerically solving both equations for Astar
    Astar1=vpasolve(Eqn2(1),Astar); 
    Astar2=vpasolve(Eqn2(2),Astar);
    % Compiling both results for Astar and selecting the smallest value,
    % which yields the correct solution for Astar
    Astar1=[Astar1;Astar2];
    AstarP=min(Astar1);
    K=find(Astar1==min(Astar1));
    Sstar=subs(Sstar,Astar,AstarP);
    Astar=AstarP;
    % Calculating fluorescence rates
    FL=kfl*((0.25*ktt*(Astar.^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))./(kfret*(S0-Sstar)+kfl+kNR));
    Aastar=((0.25*ktt*(Astar.^2)+krisc*0.75*(ktt/(kic+krisc))*(Astar^2))./(kfret*(S0-Sstar)+kfl+kNR));
    % The matrices below store all relavant emission data that is to be
    % plotted
    FLA=[FLA,FL(1)*(kext/(kext+kexP*(S0-Sstar(K))))];
    PH=[PH,kph*Sstar(K)];
    EmR=[EmR,(kph*Sstar(K))./(FL(1)*(kext/(kext+kexP*(S0-Sstar(K)))))];
    Eff=[Eff,(FL(1)*(kext/(kext+kexP*(S0-Sstar(K)))))/(kex*I*(S0-Sstar(K)))];
    EffEx=[EffEx,FL(1)/I];
end