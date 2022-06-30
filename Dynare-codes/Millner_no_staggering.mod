/*
*
* 
Last modification : 27/06/2022
Author : Davi Méaille 
*
*
Description : 
*
This mod.file implements the Dietz and Millner model in a RBC model with deviation in trend from 
steady-state temperature. It models the impact of a shock of temperature and the 
response in building of public adaptation capital at different stages of climate change. 
We assume a time-to-build hypothesis on the building of public capital, as in Bouakez and Guillard (2020). 
*
*
We use symbolic functions to compute the steady-state for protection capital and the following 
steady-state for output. 
*
This mod_file allows to implement all of the different specifications steps by steps with the different 
parameters at the beginning. 
*
*
*
*/ 


/* 
    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            I - INTRODUCTION AND DECLARATION OF VARIABLES
    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/

//------// Utility function form //-----//
// Set 1 for log additive utility (log of multiplicative), 2 for multiplicative utility, or some other values for additive CRRA utility 
@#define log_utility = 2

//-----// Step by step implementation //-----//
// Set 10 for protective without adaptation capital, 11 for both, 00 for none, 01 for only adaptation
@#define productive_protection = 11
@#define productive_capital = 1
@#define protection_capital = 01

// Set 1 for no time to build, 4 for one year (4 quarters) and 12 for three years (12 quarters)
@#define time_to_build = 1

// Share of adaptative capital kP over Y at the steady-state
share = 0.30;


var 

// General variables
Y  ${Y}$           (long_name='GDP')
C  ${C}$           (long_name='Consumption')
n  ${N}$           (long_name='Hours worked')
k  ${K}$           (long_name='Private capital')
Invest ${I}$       (long_name='Private investment')
q  ${q}$           (long_name='Tobin q for private investment')
Uc ${U_C}$         (long_name='Marginal utility of consumption')
Un ${U_N}$         (long_name='Marginal disutility of labour')
W  ${W}$           (long_name='Wage of the household')
R  ${F_K}$           (long_name='Interest rate')
c  ${\frac{C}{Y}}$      (long_name='Steady-state ratio for private consumption')
in ${\frac{I}{Y}}$      (long_name='Steady-state ratio for private investment')


// Productive capital
    @#if productive_capital == 1
kI ${K_{G,I}}$     (long_name='Public productive capital')
AI ${A_{G,I}}$     (long_name='Public investment in productive capital')
GI ${G_{I}}$       (long_name='Public effective investment in productive capital')
qI ${q_I}$         (long_name='Tobin q for public investment in productive capital')
ai ${\frac{AI}{Y}}$(long_name='Steady-state ratio for public productive investment')
    @#endif

    @#if protection_capital == 1
// Protection capital 
kP ${K_{G,P}}$     (long_name='Public protection capital')
AP ${A_{G,P}}$     (long_name='Public investment in protection capital')
GP ${G_{P}}$       (long_name='Public effective investment in protection capital')
qP ${q_P}$         (long_name='Tobin q for public investment in protection capital')
ap ${\frac{AP}{Y}}$(long_name='Steady-state ratio for public protection investment')
    @#endif

// Weather shock
temp ${\tau^W}$    (long_name='Function of temperature damage')
damage ${D}$       (long_name='Function of damage on the economy, net from adaptation capital')

@#if productive_protection == 01 || productive_protection == 11
// Debt
T  ${T}$           (long_name='Lump-sum taxes')
B  ${B}$           (long_name='Government bond')
S  ${S}$           (long_name='Primary surplus')
Rb ${R_B}$         (long_name='Interest rate for debt')
G  ${G}$           (long_name='Public expenditures')
g  ${\frac{G}{Y}}$           (long_name='Public expenditures over GDP')

// Debt
b  ${\frac{B}{Y}}$           (long_name='Debt-output ratio')
s  ${\frac{S}{Y}}$           (long_name='Primary surplus to income ratio')
t  ${\frac{T}{Y}}$ (long_name='Ratio of lump-sum tax over output')
@#endif
; 



varexo 
epsilon ${\varepsilon_\tau}$   (long_name='Temperature shock')
;



parameters 
// General parameters
betta    ${\beta}$            (long_name='Discount factor')
khi      ${\khi}$             (long_name='Coefficient of labour disutility')
sigmaC   ${\sigma_C}$         (long_name='Inverse intertemporal elasticity of substitution')
sigmaN   ${\sigma_N}$         (long_name='Inverse Frisch elasticity of labor')
alfa     ${\alpha}$           (long_name='Capital share')
alfaG    ${\alpha_G}$         (long_name='Capital share of public productive capital')
deltta   ${\delta}$           (long_name='Depreciation rate of private capital')
delttaI  ${\delta_I}$         (long_name='Depreciation rate of public productive capital')
delttaP  ${\delta_P}$         (long_name='Depreciation rate of public protection capital')
omga     ${\omega}$           (long_name='Coefficient of investment friction for household')
psi      ${\psi}$             (long_name='Coefficient of investment friction for government capital')
N 
gama     ${\gamma}$           (long_name='Preference parameter')
sigma    ${\sigma}$           (long_name='Other preference parameter') 

// Steady-state values 
Yss         ${Y_{SS}}$        (long_name='Initial value for the steady-state of GDP')
css         ${C_{SS}}$        (long_name='Initial value for the steady-state of households consumption')
nss         ${N_{SS}}$        (long_name='Initial value for the steady-state of households work over leisure ratio')
kss         ${K_{SS}}$        (long_name='Initial value for the steady-state of households capital')
kIss        ${K_{G,I,SS}}$    (long_name='Initial value for the steady-state of government productive capital')
Investss    ${I_{SS}}$        (long_name='Initial value for the steady-state of households investment') 
AIss        ${A_{I,SS}}$      (long_name='Initial value for the steady-state of government productive capital')
GIss        ${G_{I,SS}}$      (long_name='Initial value for the steady-state of effective government productive investment') 
qss         ${q_{SS}}$        (long_name='Initial value for the steady-state of the Tobin q for private capital')
qIss        ${q_{I,SS}}$      (long_name='Initial value for the steady-state of the Tobin q for government productive investment')
qPss        ${q_{P,SS}}$      (long_name='Initial value for the steady-state of the Tobin q for government protection investment') 
kPss        ${K_{G,P,SS}}$    (long_name='Initial value for the steady-state of government protection capital')
GPss        ${G_{P,SS}}$      (long_name='Initial value for the steady-state of effective government protection investment')
APss        ${A_{P,SS}}$      (long_name='Initial value for the steady-state of government protection investment')
damagess    ${\Omega_{SS}}$   (long_name='Initial value for the steady-state of temperature damage')
tempss      ${\tau^W_{SS}}$   (long_name='Steady-state value of temperature')

// Debt parameters
Tss         ${T_{SS}}$        (long_name='Steady-state value for lump-sum taxes') 
bss         ${b_{SS}}$        (long_name='Steadu-state value for the ratio of debt over output')
sss         ${s_{SS}}$        (long_name='Long-run value of debt over GDP ratio')
lambda      ${\lambda}$       (long_name='Reaction of the primary surplus to income ratio')
mu          ${\mu}$           (long_name='Intercept of the rule of primary surplus')

// Protection capital parameters
b1  ${b_1}$         (long_name='Gross damage multiplier parameter')
b2  ${b_2}$         (long_name='Gross damage multiplier parameter')        
a1  ${a_1}$         (long_name='Residual damage multiplier parameter - effectiveness of adaptation')
a2  ${a_2}$         (long_name='Residual damage multiplier parameter - effectiveness of adaptation')
a3  ${a_3}$         (long_name='Residual damage multiplier parameter - effectiveness of adaptation')

// Weather parameters
tempopt  ${\tau^*}$     (long_name='Optimal temperature')
tempm    ${\tau_M}$     (long_name='Mean of actual temperature')
rhot     ${\rho_\tau}$  (long_name='Auto-correlation of temperature shock')
;




/* 
        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                II - Parametrization
        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/


// friction parameters
omga = 2;
psi = 2;


// general parameters
alfaG = 0.1;
alfa = 0.4;
betta = 0.98;
deltta = 0.06;
delttaI = 0.06;
delttaP = 0.02;
sigmaC = 1.5;
sigmaN = .5;
rhot = 0.95;
sigma = 2;


// damage parameters
b1 = 0.115;
a1 = 0.0012; 
a2 = 0.023;
a3 = 2;
tempm = 36;
tempopt = 33; 



// Steady state parameters depending on if we take public capital and protection capital into account or not

qss = 1; 
qPss = 1;
qIss = 1;
tempss = tempm;
nss = 1/3;

@#if time_to_build == 1
N = 1;
@#elseif time_to_build == 4
N = 4;
@#else 
N = 12;
@#endif


    @#if productive_protection == 00 
damagess = 1/(1 + a1 * (tempss - tempopt) + a2 * (tempss - tempopt)^2);
Yss =  (  (damagess * nss)^alfa * (betta * (1 - alfa) / (qss * (1 - betta + betta * deltta)))^(1 - alfa) )^(1 / alfa);
kss = betta * (1 - alfa) * Yss / (qss * (1 - betta + betta * deltta)) ; 
Investss = deltta * kss;
css = Yss - Investss; 

syms gam;
eqn = gam/(1 - gam) == css/Yss * nss/(alfa*(1 - nss));
val1 = vpasolve(eqn, gam);
gama = val1;

    @#endif

    @#if productive_protection == 01
%------------------------------------------------------------------------------------------------------------
% Computation of the value of the steady-state for ${Y}$ and ${K_{G,P}}$ with matlab symbolic equation solver
%------------------------------------------------------------------------------------------------------------
syms Y C k Invest kP AP gam b damage;
assume(Y > 0 & kP > 0 & b > 0 );

eqn1 = kP/Y == share;
eqn2 = Y ==  ((1 + b1 * (kP)^b)/(1 + b1 * (kP)^b + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2) * nss)^alfa * (k)^(1 - alfa)  ; 
eqn3 = Y == C + Invest + AP;
eqn4 = gam/(1 - gam) == C/Y * nss/(alfa*(1 - nss));
% Private capital 
eqn5 = Invest == deltta * k; 
eqn6 = k == betta * (1 - alfa) * Y / (qss * (1 - betta * (1 - deltta))); 
eqn9 = AP == delttaP * kP;
eqn10 = ...
betta^N * alfa * (b1 * b * kP^(b - 1) * (a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2)/ ( 1 + b1 * kP^(b) + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2 )^2)* damage^(alfa - 1) ...
    * ( (nss)^(alfa) * k^(1 - alfa) ) ...
    == qPss - qPss * (betta * (1 - delttaP));
eqn11 = damage == (1 + b1 * (kP)^b)/(1 + b1 * (kP)^b + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2);


[val1, val2, val3, val4, val5, val8, val9, val10, val11] = vpasolve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn9, eqn10, eqn11], ...
                   [b,  Y, C, k, Invest, kP, AP,  gam, damage], ...
                   [0.5,11 ,7,70, 5,  1.5, 0.09,  0.43, 0.98]);

b2 = val1;
Yss = val2;
css = val3;
kss = val4;
Investss = val5;
kPss = val8;
APss = val9;
gama = val10;
damagess = val11;


apss = APss / Yss;
    @#endif



    @#if productive_protection == 10
damagess = 1/(1 + a1 * (tempss - tempopt) + a2 * (tempss - tempopt)^2);
Yss = ( damagess * (nss)^alfa * (betta * (1 - alfa) / (qss * (1 - betta + betta * deltta)))^(1 - alfa)  * (betta^(N) * (alfaG) / (qIss * (1 - betta + betta * delttaI)))^(alfaG) )^(1 / (alfa - alfaG));
kIss = betta^(N) * alfaG * Yss / (qIss * (1 - betta + betta * delttaI));
AIss = delttaI * kIss;
GIss = AIss;
kss = betta * (1 - alfa) * Yss / (qss * (1 - betta + betta * deltta)) ; 
Investss = deltta * kss;
css = Yss - Investss - GIss;
    @#endif


    @#if productive_protection == 11

syms Y C k Invest kI AI kP AP gam b damage;
assume(Y > 0 & kP > 0 & b > 0 );

% Calibration de b2
eqn1 = kP/(kI + kP) == share;
% General equations
eqn2 = Y ==  ((1 + b1 * (kP)^b)/(1 + b1 * (kP)^b + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2) * nss)^alfa * (k)^(1 - alfa)  * (kI)^(alfaG) ; 
eqn3 = Y == C + Invest + AP + AI;
eqn4 = gam/(1 - gam) == C/Y * nss/(alfa*(1 - nss));
% Private capital 
eqn5 = Invest == deltta * k; 
eqn6 = k == betta * (1 - alfa) * Y / (qss * (1 - betta * (1 - deltta))); 
% Public capital 
eqn7 = AI == delttaI * kI;
eqn8 = kI == betta^(N) * alfaG * Y / (qIss * (1 - betta + betta * delttaI));
% Protection capital 
eqn9 = AP == delttaP * kP;
eqn10 = ...
betta^N * alfa * (b1 * b * kP^(b - 1) * (a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2)/ ( 1 + b1 * kP^(b) + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2 )^2)* damage^(alfa - 1) ...
    * ( (nss)^(alfa) * k^(1 - alfa) * kI^(alfaG) ) ...
    == qPss - qPss * (betta * (1 - delttaP));
eqn11 = damage == (1 + b1 * (kP)^b)/(1 + b1 * (kP)^b + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2);



[val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11] = vpasolve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10, eqn11], ...
                   [b,  Y, C, k, Invest, kI, AI, kP, AP,  gam, damage], ...
                   [0.5,11 ,7,70, 5, 24, 0.48,  1.5, 0.09,  0.43, 0.98]);

b2 = val1;
Yss = val2;
css = val3;
kss = val4;
Investss = val5;
kIss = val6;
AIss = val7;
kPss = val8;
APss = val9;
gama = val10;
damagess = val11;


apss = APss / Yss;
    @#endif


// For all utility cases

khi = alfa * Yss / (css^(1/sigmaC) * nss^((1+sigmaN)/sigmaN));




// For debt

Rss = (1 - alfa) * Yss/kss;
Rbss = 1/betta - 1;
sss = 0.6 * Rbss;

@#if productive_protection == 11
Tss = AIss + APss + sss * Yss;
bss = sss/Rbss;
lambda = 0.05;
mu = (Yss * (Rbss - lambda) * 0.6 + 2 * (APss + AIss)) / Yss;
@#elseif productive_protection == 01
Tss = APss + sss * Yss;
bss = sss/Rbss;
lambda = 0.05;
mu = (Yss * (Rbss - lambda) * 0.6 + 2 * (APss)) / Yss;
@#endif 

khi = alfa * Yss / (css^(1/sigmaC) * nss^((1+sigmaN)/sigmaN));




/* 
        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                III - DECLARATION OF THE MODEL
        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/

model; 
%%%%------- General equations --------%%%%%

[name='Interest rate']
R = (1 - alfa) * Y/k(-1);

[name='Wage']
W = alfa * Y/n;

[name='Steady-state ratio for private consumption']
c = C/Y;

    @#if log_utility == 1  
[name='Marginal utility of consumption']
Uc = gama/C;

[name='Marginal disutility of labour']
Un = -(1 - gama)/(1 - n);

    @#elseif log_utility == 2
[name='Marginal utility of consumption']
Uc = gama * C^(gama * (1 - sigma) - 1) * (1 - n)^((1 - gama)*(1 - sigma));

[name='Marginal disutility of labour']
Un = - (1 - gama) * C^(gama * (1 - sigma)) * (1 - n)^((1 - gama)*(1 - sigma) - 1);

    @#else
[name='Marginal utility of consumption']
Uc = (C)^(-1/sigmaC);

[name='Marginal disutility of labour']
Un = - khi * n^(1/sigmaN);

    @#endif




    @#if productive_capital == 1
[name='Production of the firm']
Y = (damage * n)^(alfa) * k(-1)^(1 - alfa) * kI(-1)^(alfaG); 
    @#else
Y = (damage * n)^(alfa) * k(-1)^(1 - alfa);
    @#endif

[name='FCO working hours']
- Un/Uc = W; 




%%%%-------- Private capital ----------%%%%%
[name='Capital motion']
k = (1-deltta)*k(-1)+ (1-(omga/2*(Invest/Invest(-1)-1)^2))*Invest;   

[name='FOC in private capital']
betta * (1 - alfa) * Y(+1) / k =
    q * Uc/Uc(+1) - q(+1) * betta * (1 - deltta); 

[name='FCO on investment'] 
q * (1 - omga/2*((Invest/Invest(-1))-1)^2 - omga*((Invest/Invest(-1))-1)*Invest/Invest(-1))
 +betta*Uc(+1)/Uc*q(+1)*omga*((Invest(+1)/Invest)-1)*(Invest(+1)/Invest)^2 = 1;

[name='Steady-state ratio for private investment']
in = Invest/Y;



%%%%-------- Resource constraint of the economy -------%%%%

    @#if productive_protection == 00
[name='Resource constraint of the economy']
Y = C + Invest;
    @#endif

    @#if productive_protection == 01 
[name='Resource constraint of the economy']
Y = C + Invest + GP;
    @#endif

    @#if productive_protection == 10
[name='Resource constraint of the economy']
Y = C + Invest + GI;
    @#endif

    @#if productive_protection == 11
[name='Resource constraint of the economy']
Y = C + Invest + GP + GI;
    @#endif







%%%%%------------ Equations related to debt -----------%%%%%

    @#if productive_protection == 11
[name='Primary surplus to output ratio']
s = (T - GI - GP)/Y ;

[name='Primary surplus']
S = T - GI - GP;

[name='Government budget constraint']
B(-1) * (1 + Rb(-1)) = s * Y + B;

[name='Adjustment of fiscal policy from Leeper(2010)']
T/Y = lambda * b(-1) + mu - (GI + GP)/Y;

[name='Debt-output ratio']
b = B/Y;

[name='Interest rate for debt']
1 = betta * Uc(+1)/Uc * (1 + Rb);

[name='Ratio of lump-sum tax over GDP']
t = T/Y;



    @#elseif productive_protection == 01

[name='Primary surplus to output ratio']
s = (T - GP)/Y ;

[name='Primary surplus']
S = T - GP;

[name='Government budget constraint']
B(-1) * (1 + Rb(-1)) = s * Y + B;

[name='Adjustment of fiscal policy from Leeper(2010)']
T/Y = lambda * b(-1) + mu - ( GP)/Y;

[name='Debt-output ratio']
b = B/Y;

[name='Interest rate for debt']
1 = betta * Uc(+1)/Uc * (1 + Rb);

[name='Ratio of lump-sum tax over GDP']
t = T/Y;

    @#endif



%%%%-------- Public capital ----------%%%%%
    @#if productive_capital == 1

[name='FCO on productive investment AI']
qI * (1 -  psi/2*((AI/AI(-1))-1)^2 -  psi*((AI/AI(-1))-1)*AI/AI(-1))
    +betta*Uc(+1)/Uc*qI(+1)* psi*((AI(+1)/AI)-1)*(AI(+1)/AI)^2 = 1;

[name='Steady-state ratio for public productive investment']
ai = AI/Y;

    @#endif 

    @#if productive_capital == 1
        @#if time_to_build == 1
[name='Capital motion for productive capital']
kI = (1-delttaI)*kI(-1)+(1-( psi/2*(AI/AI(-1)-1)^2))*AI;

[name='FCO on public productive capital']
betta^(1) * alfaG * Y(+1) / kI =
    qI * Uc/Uc(+1) - qI(+1) * (betta * (1 - delttaI)) * Uc(+1)/Uc(+1);

        @#elseif time_to_build == 4
[name='Capital motion for productive capital']
kI = (1-delttaI)*kI(-1)+(1-( psi/2*(AI(-3)/AI(-4)-1)^2))*AI(-3);

[name='FCO on public productive capital']
betta^(4) * alfaG * Y(+4) / kI(+3) =
    qI * Uc/Uc(+4) - qI(+1) * (betta * (1 - delttaI)) * Uc(+1)/Uc(+4);

        @#else 
[name='Capital motion for productive capital']
kI = (1-delttaI)*kI(-1)+(1-( psi/2*(AI(-11)/AI(-12)-1)^2))*AI(-11);

[name='FCO on public productive capital']
betta^(12) * alfaG * Y(+12) / kI(+11) =
    qI * Uc/Uc(+12) - qI(+1) * (betta * (1 - delttaI)) * Uc(+1)/Uc(+12);
    @#endif


[name='Public investment in productive capital']              
GI = AI;

  @#endif


%%%%-------- Protection capital ----------%%%%%

    @#if protection_capital == 1

[name='FCO on protection investment'] 
qP * (1 -  psi/2*((AP/AP(-1))-1)^2 -  psi*((AP/AP(-1))-1)*AP/AP(-1))
    + betta*Uc(+1)/Uc*qP(+1)* psi*((AP(+1)/AP)-1)*(AP(+1)/AP)^2 = 1;

[name='Steady-state ratio for public protection investment']
ap = AP/Y;

    @#endif


// Capital motion

    @#if protection_capital == 1

        @#if time_to_build == 1
[name='Capital motion for protection capital']
kP = (1-delttaP)*kP(-1)+(1-( psi/2*(AP/AP(-1)-1)^2))*AP;

        @#elseif time_to_build == 4
[name='Capital motion for protection capital']
kP = (1-delttaP)*kP(-1)+(1-( psi/2*(AP(-3)/AP(-4)-1)^2))*AP(-3);

        @#else
[name='Capital motion for protection capital']
kP = (1-delttaP)*kP(-1)+(1-( psi/2*(AP(-11)/AP(-12)-1)^2))*AP(-11);

    @#endif

[name='Public investment in protection capital']
GP = AP;

@#endif

%-----------------------------------------------------

// FOC

// Protection capital without productive capital

    @#if productive_protection == 01

        @#if time_to_build == 1
[name='FCO on protection capital']
betta^1 * b1 * alfa * b2 * kP^(b2 - 1) * (a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2) * damage^(alfa - 1) * ( (n(+1))^(alfa) * k^(1 - alfa) ) 
    / ( 1 + b1 * kP^(b2) + a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2 )^2  
    = qP * Uc/Uc(+1) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+1);

        @#elseif time_to_build == 4
[name='FCO on protection capital']
betta^4 * b1 * alfa * b2 * kP(3)^(b2 - 1) * (a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2) * damage^(alfa - 1) *  ( (n(+4))^(alfa) * k(3)^(1 - alfa) ) 
    / ( 1 + b1 * kP(3)^(b2) + a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2 )^2  
    = qP * Uc/Uc(+4) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+4);

        @#else

[name='FCO on protection capital']
betta^12 * alfa * b1 * b2 * kP(11)^(b2 - 1) * (a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2) * damage^(alfa - 1) * ( (n(+12))^(alfa) * k(11)^(1 - alfa) ) 
    / ( 1 + b1 * kP(11)^(b2) + a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2 )^2  
    = qP * Uc/Uc(+12) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+12);

    @#endif

  @#endif

// Protection capital with productive capital 

    @#if productive_protection == 11

        @#if time_to_build == 1
[name='FCO on protection capital']
betta^1 * alfa * (b1 * b2 * kP^(b2 - 1) * (a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2)/ ( 1 + b1 * kP^(b2) + a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2 )^2 ) 
    * damage^(alfa - 1) 
    * ( (n(+1))^(alfa) * k^(1 - alfa) * kI^(alfaG) )  
    = qP * Uc/Uc(+1) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+1);

        @#elseif time_to_build == 4

[name='FCO on protection capital']
betta^4 * alfa * (b1 * b2 * kP(3)^(b2 - 1) * (a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2) / ( 1 + b1 * kP(3)^(b2) + a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2 )^2 )
    *damage^(alfa - 1) 
    * ( (n(+4))^(alfa) * k(3)^(1 - alfa) * kI(3)^(alfaG) )  
    = qP * Uc/Uc(+4) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+4);

        @#else
[name='FCO on protection capital']
betta^12 * alfa * (b1 * b2 * kP(11)^(b2 - 1) * (a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2)/ ( 1 + b1 * kP(11)^b2 + a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2 )^2)
    * damage^(alfa - 1) 
    * ( (n(+12))^(alfa) * k(11)^(1 - alfa) * kI(11)^(alfaG))   
    = qP * Uc/Uc(+12) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+12);

    @#endif
  @#endif


%---------------------------------------------------






%%%% ------- Natural disaster --------%%%%

    @#if protection_capital == 1
[name='Damage function with protection capital']
damage = (1 + b1 * (kP(-1))^b2)/(1 + b1 * (kP(-1))^b2 + a1 * (temp - tempopt) + a2 * (temp - tempopt)^2);
    @#else
[name='Damage function']
damage = 1/(1 + a1 * (temp - tempopt) + a2 * (temp - tempopt)^2);
    @#endif

[name='Climatic shock']
temp = (1 - rhot) * tempm + rhot * temp(-1)  + epsilon;


@#if productive_protection == 01
G = GP;
g = G/Y;
@#elseif productive_protection == 11
G = GI + GP;
g = G/Y;
@#endif



end; 


/* 
        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                IV - DECLARATION OF THE STEADY STATE
        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/


steady_state_model; 

n = nss;
Y = Yss;
k = kss; 
Invest = deltta * kss;
q = qss;
C = css;
temp = tempm;
damage = damagess;
W = alfa * Y/n;
R = (1 - alfa) * Y/k;
c = C/Y;
in = Invest/Y;


/*          IMPLEMENTATION         */
    @#if productive_capital == 1
kI = kIss; 
AI = delttaI * kIss;
GI = AI;
qI = qIss;
ai = AI/Y;
    @#endif 

    @#if protection_capital == 1
kP = kPss;
AP = delttaP * kPss; 
GP = AP;
qP = qPss;
ap = AP/Y;

    @#endif




    @#if log_utility == 1
Uc = gama/C;
Un = - (1 - gama)/(1 - n);
    @#elseif log_utility == 2
Uc = gama * C^(gama * (1 - sigma) - 1) * (1 - n)^((1 - gama)*(1 - sigma));
Un = - (1 - gama) * C^(gama * (1 - sigma)) * (1 - n)^((1 - gama)*(1 - sigma) - 1);
    @#else
Uc = (C)^(-1/sigmaC);
Un = - khi * n^(1/sigmaN);
    @#endif




/*         DEBT          */
        @#if productive_protection == 11
Rb = 1/betta - 1;
B = (mu * Yss - 2 * GI - 2 * GP)/(Rb - lambda);
b = B/Y;
T = lambda * B + mu * Yss - GI - GP;
s = (T - GI - GP)/Y;
G = GI + GP;
t = T/Y;
S = s * Y;
g = G/Y;
        @#elseif productive_protection == 01
Rb = 1/betta - 1;
B = (mu * Yss - 2 * GP)/(Rb - lambda);
b = B/Y;
T = lambda * B + mu * Yss - GP;
s = (T - GP)/Y;
G = GP;
t = T/Y;
S = s * Y;
g = G/Y;
        @#endif 

end; 


/* 
                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        V - SIMULATIONS
                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/


shocks;
var epsilon; stderr .01;
end;

write_latex_dynamic_model;

resid(2);
steady; 
check;
options_.TeX=1;
@#if productive_protection == 01
stoch_simul(nomoments, nocorr, nofunctions, irf = 20) Y C n k Invest W R kP AP G B T S Rb  c in ap g b t s;
@#elseif productive_protection == 10
stoch_simul(nomoments, nocorr, nofunctions, irf = 100) Y c n k in W R kI ai;
@#elseif productive_protection == 11 
stoch_simul(nomoments, nocorr, nofunctions, irf = 120) Y C n k Invest W R kI AI kP AP G B T S Rb  c in ai ap g b t s;
@#else 
stoch_simul(nomoments, nocorr, nofunctions, irf = 50) Y C n k Invest W R c in q;
@#endif

names = M_.endo_names;
tex = M_.endo_names_tex;

