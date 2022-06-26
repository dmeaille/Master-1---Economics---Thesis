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
This allows to draw different impulse-response functions, that are generally smoother that in the Leeper(2010) case. 
*
There is only one recquirement : for the Bouakez and Guillard(2020) time-to-build specification to work, 
when included, we have to add friction on investment capital, for private, public productive and 
protection capital. 
*
This allows also to check for the robustness of our paper's results, as the time-to-build is an important 
explaining factor.  
*
This .mod file implements also the calibration for the parameter b2 to calibrate the share of adaptative 
capital over the total of public investment. 
*
Eventually, this .mod file implements different response fonction for public capital, from the optimal 
response in investment to other ad hoc policies.
*
*
*
*/ 


/* 
    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            I - INTRODUCTION AND DECLARATION OF VARIABLES
    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/

// Parametrization
// Set 1 for log additive utility (log of multiplicative), 2 for multiplicative utility, or some other values for additive CRRA utility 
@#define log_utility = 1

// Set 1 for no time to build, 4 for one year (4 quarters) and 12 for three years (12 quarters)
@#define time_to_build = 1

// Share of adaptative capital kP over total public capital kP + kI
share = 0.15;

// Set 1 for optimal response of adaptative capital, 2 for no response (constant stock) and 3 for constant share of investment over PIB. 
@#define optimal = 3

// Set 1 for heatwave, 2 for flood, 0 for a shock affecting the total factor productivity
@#define type = 1

var 
        /*§§§§§§§§§ Variable en niveau §§§§§§§§§§*/
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
R  ${R}$           (long_name='Interest rate')
G  ${G}$           (long_name='Public expenditures')

// Debt
T  ${T}$           (long_name='Lump-sum taxes')
B  ${B}$           (long_name='Government bond')
S  ${S}$           (long_name='Primary surplus')
Rb ${R_b}$         (long_name='Interest rate for debt')

// Productive capital
kI ${K_{G,I}}$     (long_name='Public productive capital')
AI ${A_{G,I}}$     (long_name='Public investment in productive capital')
GI ${G_{I}}$       (long_name='Public effective investment in productive capital')
qI ${q_I}$         (long_name='Tobin q for public investment in productive capital')

// Protection capital 
kP ${K_{G,P}}$     (long_name='Public protection capital')
AP ${A_{G,P}}$     (long_name='Public investment in protection capital')
GP ${G_{P}}$       (long_name='Public effective investment in protection capital')
qP ${q_P}$         (long_name='Tobin q for public investment in protection capital')
kp ${\frac{K_P}{Y}}$(long_name='Steady-state ratio for public adaptation investment')

// Weather shock
damage ${D}$       (long_name='Function of damage on the economy, net from adaptation capital')
share   ${Share}$  (long_name='Share of adaptative capital over total public capital')

 
temp ${\tau^W}$    (long_name='Function of temperature damage')


        /*$$$$$$$$$ Variable en ratio GDP §§§§§§§§§*/


// General variables
c  ${\frac{C}{Y}}$      (long_name='Steady-state ratio for private consumption')
in ${\frac{I}{Y}}$      (long_name='Steady-state ratio for private investment')
g  ${g}$                (long_name='Public expenditures over GDP')

// Debt
t  ${t}$           (long_name='Lump-sum taxes ratio over GDP')
b  ${b}$           (long_name='Debt-output ratio')
s  ${s}$           (long_name='Primary surplus to income ratio')
t  ${\frac{T}{Y}}$ (long_name='Ratio of lump-sum tax over output')

// Public capital
ai ${\frac{A_I}{Y}}$    (long_name='Steady-state ratio for public productive investment')
ap ${\frac{A_P}{Y}}$    (long_name='Steady-state ratio for public protection investment')

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
Yss         ${Y_{SS}}$        (long_name='Steady-state value for the steady-state of GDP')
css         ${C_{SS}}$        (long_name='Steady-state value for the steady-state of households consumption')
nss         ${N_{SS}}$        (long_name='Steady-state value for the steady-state of households work over leisure ratio')
kss         ${K_{SS}}$        (long_name='Steady-state value for the steady-state of households capital')
kIss        ${K_{G,I,SS}}$    (long_name='Steady-state value for the steady-state of government productive capital')
Investss    ${I_{SS}}$        (long_name='Steady-state value for the steady-state of households investment') 
AIss        ${A_{I,SS}}$      (long_name='Steady-state value for the steady-state of government productive capital')
GIss        ${G_{I,SS}}$      (long_name='Steady-state value for the steady-state of effective government productive investment') 
qss         ${q_{SS}}$        (long_name='Steady-state value for the steady-state of the Tobin q for private capital')
qIss        ${q_{I,SS}}$      (long_name='Steady-state value for the steady-state of the Tobin q for government productive investment')
qPss        ${q_{P,SS}}$      (long_name='Steady-state value for the steady-state of the Tobin q for government protection investment') 
kPss        ${K_{G,P,SS}}$    (long_name='Steady-state value for the steady-state of government protection capital')
GPss        ${G_{P,SS}}$      (long_name='Steady-state value for the steady-state of effective government protection investment')
APss        ${A_{P,SS}}$      (long_name='Steady-state value for the steady-state of government protection investment')
damagess    ${\Omega_{SS}}$   (long_name='Steady-state value for the steady-state of temperature damage')
tempss      ${\tau^W_{SS}}$   (long_name='Steady-state value of temperature')
apss        ${\frac{AP_{SS}}{Y_{SS}}$       (long_name='Steady-state ratio of adaptative capital')
Rss         ${R_{SS}}$        (long_name='Steady-state value for the real interest rate')

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
sigmaC = 1.75;
sigmaN = 2;
rhot = 0.95;
sigma = 2;


// damage parameters
b1 = 0.115;
a1 = 0.0012; 
a2 = 0.023;
a3 = 2;
tempm = 30;
tempm2 = 27; 
tempopt = 27;






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



syms Y C k Invest kI AI kP AP gam b damage;
assume(Y > 0 & kP > 0 & b > 0 );

% Calibration de b2
eqn1 = kP/(kI + kP) == share;
% General equations

    @#if type == 0
eqn2 = Y == (1 + b1 * (kP)^b)/(1 + b1 * (kP)^b + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2) * (nss)^alfa * (k)^(1 - alfa)  * (kI)^(alfaG) ; 
eqn3 = Y == C + Invest + AP + AI;
eqn4 = gam/(1 - gam) == C/Y * nss/(alfa*(1 - nss));
% Private capital 
eqn5 = Invest == deltta * k; 
eqn6 = k == betta * (1 - alfa) * Y / (qss * (1 - betta * (1 - deltta))); 
    @#elseif type == 1
eqn2 = Y ==  ((1 + b1 * (kP)^b)/(1 + b1 * (kP)^b + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2) * nss)^alfa * (k)^(1 - alfa)  * (kI)^(alfaG) ; 
eqn3 = Y == C + Invest + AP + AI;
eqn4 = gam/(1 - gam) == C/Y * nss/(alfa*(1 - nss));
% Private capital 
eqn5 = Invest == deltta * k; 
eqn6 = k == betta * (1 - alfa) * Y / (qss * (1 - betta * (1 - deltta))); 
    @#elseif type == 2
eqn2 = Y == (nss)^alfa * ((1 + b1 * (kP)^b)/(1 + b1 * (kP)^b + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2) * k)^(1 - alfa)  * (kI)^(alfaG) ; 
eqn3 = Y == C + Invest + AP + AI;
eqn4 = gam/(1 - gam) == C/Y * nss/(alfa*(1 - nss));
% Private capital 
eqn5 = Invest == (1 - (1 - deltta) * damage) * k; 
eqn6 = k == betta * (1 - alfa) * Y / (qss * (1 - betta * (1 - deltta) * damage)); 
    @#endif

% Public capital 
eqn7 = AI == delttaI * kI;
eqn8 = kI == betta^(N) * alfaG * Y / (qIss * (1 - betta + betta * delttaI));
% Protection capital 
eqn9 = AP == delttaP * kP;

    @#if type == 0
eqn10 = ...
betta^N * b1 * b * kP^(b - 1) * (a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2) * ( (nss)^(alfa) * k^(1 - alfa) * kI^(alfaG) ) ...
    / ( 1 + b1 * kP^(b) + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2 )^2  ...
    == qPss - qPss * (betta * (1 - delttaP));
    @#elseif type == 1
eqn10 = ...
betta^N * alfa * (b1 * b * kP^(b - 1) * (a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2)/ ( 1 + b1 * kP^(b) + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2 )^2)* damage^(alfa - 1) ...
    * ( (nss)^(alfa) * k^(1 - alfa) * kI^(alfaG) ) ...
    == qPss - qPss * (betta * (1 - delttaP));
    @#elseif type == 2
eqn10 = ...
betta^N * (1 - alfa) * (b1 * b * kP^(b - 1) * (a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2)/ ( 1 + b1 * kP^(b) + a1 * (tempm - tempopt) + a2 * (tempm - tempopt)^2 )^2) * damage^(-alfa) ...
    * ( (nss)^(alfa) * k^(1 - alfa) * kI^(alfaG) ) ...
    == qPss - qPss * (betta * (1 - delttaP));
    @#endif 

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

Rss = (1 - alfa) * Yss/kss;
Rbss = 1/betta - 1;
sss = 0.6 * Rbss;
Tss = AIss + APss + sss * Yss;
apss = APss / Yss;
bss = sss/Rbss;


// For all utility cases
khi = alfa * Yss / (css^(1/sigmaC) * nss^((1+sigmaN)/sigmaN));



lambda = 0.05;
mu = (Rbss - lambda) * 0.6 + 2 * (APss + AIss) / Yss;

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



    @#if type == 0
[name='Production of the firm']
Y = damage * (n)^(alfa) * k(-1)^(1 - alfa) * kI(-1)^(alfaG); 
    @#elseif type == 1
[name='Production of the firm']
Y = (damage * n)^(alfa) * k(-1)^(1 - alfa) * kI(-1)^(alfaG); 
    @#elseif type == 2
[name='Production of the firm']
Y = (n)^(alfa) * ( damage * k(-1))^(1 - alfa) * kI(-1)^(alfaG); 
    @#endif


[name='FCO working hours']
- Un/Uc = W; 




%%%%-------- Private capital ----------%%%%%
    @#if type == 0 || type == 1
[name='Capital motion']
k = (1-deltta)*k(-1)+ (1-(omga/2*(Invest/Invest(-1)-1)^2))*Invest;

[name='FOC in private capital']
betta * (1 - alfa) * Y(+1) / k =
    q * Uc/Uc(+1) - q(+1) * betta * (1 - deltta); 

[name='FCO on investment'] 
q * (1 - omga/2*((Invest/Invest(-1))-1)^2 - omga*((Invest/Invest(-1))-1)*Invest/Invest(-1))
 +betta*Uc(+1)/Uc*q(+1)*omga*((Invest(+1)/Invest)-1)*(Invest(+1)/Invest)^2 = 1;
 
    @#elseif type == 2
[name='Capital motion']
k = (1-deltta)*damage*k(-1)+ (1-(omga/2*(Invest/Invest(-1)-1)^2))*Invest; 

[name='FOC in private capital']
betta * (1 - alfa) * Y(+1) / k =
    q * Uc/Uc(+1) - q(+1) * betta * (1 - deltta) * damage; 

[name='FCO on investment'] 
q * (1 - omga/2*((Invest/Invest(-1))-1)^2 - omga*((Invest/Invest(-1))-1)*Invest/Invest(-1))
 +betta*Uc(+1)/Uc*q(+1)*omga*((Invest(+1)/Invest)-1)*(Invest(+1)/Invest)^2 = 1;

    @#endif 




[name='Steady-state ratio for private investment']
in = Invest/Y;



%%%%-------- Resource constraint of the economy -------%%%%


[name='Resource constraint of the economy']
Y = C + Invest + GP + GI;





%%%%%------------ Equations related to debt -----------%%%%%


[name='Primary surplus to output ratio']
s = (T - GI - GP)/Y ;

[name='Primary surplus']
S = T - GI - GP;

[name='Government budget constraint']
B(-1) * (1 + Rb(-1)) = s * Y + B;

[name='Adjustment of fiscal policy from Leeper(2010)']
t = lambda * b(-1) + mu - (GI + GP)/Y;

[name='Debt-output ratio']
b = B/Y;

[name='Interest rate for debt']
1 = betta * Uc(+1)/Uc * (1 + Rb);

[name='Ratio of lump-sum tax over GDP']
t = T/Y;



%%%%-------- Public capital ----------%%%%%

[name='Public expenditures']
G = GI + GP;

[name='Public expenditures over GDP']
g = G/Y;

[name='FCO on productive investment AI']
qI * (1 -  psi/2*((AI/AI(-1))-1)^2 -  psi*((AI/AI(-1))-1)*AI/AI(-1))
    +betta*Uc(+1)/Uc*qI(+1)* psi*((AI(+1)/AI)-1)*(AI(+1)/AI)^2 = 1;

[name='Steady-state ratio for public productive investment']
ai = AI/Y;


        @#if time_to_build == 1
[name='Capital motion for productive capital']
kI = (1-delttaI)*kI(-1)+(1-( psi/2*(AI/AI(-1)-1)^2))*AI;

[name='Public investment in productive capital']              
GI = AI;

[name='FCO on public productive capital']
betta^(1) * alfaG * Y(+1) / kI =
    qI * Uc/Uc(+1) - qI(+1) * (betta * (1 - delttaI)) * Uc(+1)/Uc(+1);

        @#elseif time_to_build == 4
[name='Capital motion for productive capital']
kI = (1-delttaI)*kI(-1)+(1-( psi/2*(AI(-3)/AI(-4)-1)^2))*AI(-3);

[name='Public investment in productive capital']              
GI = AI(-3);

[name='FCO on public productive capital']
betta^(4) * alfaG * Y(+4) / kI(+3) =
    qI * Uc/Uc(+4) - qI(+1) * (betta * (1 - delttaI)) * Uc(+1)/Uc(+4);

        @#else 
[name='Capital motion for productive capital']
kI = (1-delttaI)*kI(-1)+(1-( psi/2*(AI(-11)/AI(-12)-1)^2))*AI(-11);

[name='Public investment in productive capital']              
GI = AI(-11);

[name='FCO on public productive capital']
betta^(12) * alfaG * Y(+12) / kI(+11) =
    qI * Uc/Uc(+12) - qI(+1) * (betta * (1 - delttaI)) * Uc(+1)/Uc(+12);
    @#endif









%%%%-------- Protection capital ----------%%%%%


@#if optimal == 1

[name='FCO on protection investment'] 
qP * (1 -  psi/2*((AP/AP(-1))-1)^2 -  psi*((AP/AP(-1))-1)*AP/AP(-1))
    +betta*Uc(+1)/Uc*qP(+1)* psi*((AP(+1)/AP)-1)*(AP(+1)/AP)^2 = 1;

[name='Steady-state ratio for public protection investment']
ap = AP/Y;



// Capital motion

 
        @#if time_to_build == 1
[name='Capital motion for protection capital']
kP = (1-delttaP)*kP(-1)+(1-( psi/2*(AP/AP(-1)-1)^2))*AP;

[name='Public investment in protection capital']
GP = AP;

        @#elseif time_to_build == 4
[name='Capital motion for protection capital']
kP = (1-delttaP)*kP(-1)+(1-( psi/2*(AP(-3)/AP(-4)-1)^2))*AP(-3);

[name='Public investment in protection capital']
GP = AP(-3);

        @#else
[name='Capital motion for protection capital']
kP = (1-delttaP)*kP(-1)+(1-( psi/2*(AP(-11)/AP(-12)-1)^2))*AP(-11);

[name='Public investment in protection capital']
GP = AP(-11);
        @#endif


// FOC
@#if type == 0
        @#if time_to_build == 1
[name='FCO on protection capital']
betta^1 * b1 * b2 * kP^(b2 - 1) * (a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2) * ( (n(+1))^(alfa) * k^(1 - alfa) * kI^(alfaG) ) 
    / ( 1 + b1 * kP^(b2) + a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2 )^2  
    = qP * Uc/Uc(+1) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+1);

        @#elseif time_to_build == 4

[name='FCO on protection capital']
betta^4 * b1 * b2 * kP(3)^(b2 - 1) * (a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2) * ( (n(+4))^(alfa) * k(3)^(1 - alfa) * kI(3)^(alfaG) ) 
    / ( 1 + b1 * kP(3)^(b2) + a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2 )^2  
    = qP * Uc/Uc(+4) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+4);

        @#else
[name='FCO on protection capital']
betta^12 * b1 * b2 * kP(11)^(b2 - 1) * (a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2) * ( (n(+12))^(alfa) * k(11)^(1 - alfa) * kI(11)^(alfaG)) 
    / ( 1 + b1 * kP(11)^b2 + a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2 )^2  
    = qP * Uc/Uc(+12) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+12);

    @#endif

@#elseif type == 1

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



@#elseif type == 2

        @#if time_to_build == 1
[name='FCO on protection capital']
betta^1 * (1 - alfa) * (b1 * b2 * kP^(b2 - 1) * (a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2) / ( 1 + b1 * kP^(b2) + a1 * (temp(+1) - tempopt) + a2 * (temp(+1) - tempopt)^2 )^2)
    *damage^(-alfa) 
    * ( (n(+1))^(alfa) * k^(1 - alfa) * kI^(alfaG) )   
    = qP * Uc/Uc(+1) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+1);

        @#elseif time_to_build == 4

[name='FCO on protection capital']
betta^4 * (1 - alfa) * (b1 * b2 * kP(3)^(b2 - 1) * (a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2) / ( 1 + b1 * kP(3)^(b2) + a1 * (temp(+4) - tempopt) + a2 * (temp(+4) - tempopt)^2 )^2)
    *damage^(-alfa) 
    * ( (n(+4))^(alfa) * k(3)^(1 - alfa) * kI(3)^(alfaG) ) 
      
    = qP * Uc/Uc(+4) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+4);

        @#else
[name='FCO on protection capital']
betta^12 * (1 - alfa) * (b1 * b2 * kP(11)^(b2 - 1) * (a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2) / ( 1 + b1 * kP(11)^b2 + a1 * (temp(+12) - tempopt) + a2 * (temp(+12) - tempopt)^2 )^2 )
    *damage^(-alfa) 
    * ( (n(+12))^(alfa) * k(11)^(1 - alfa) * kI(11)^(alfaG)) 
    = qP * Uc/Uc(+12) - qP(+1) * (betta * (1 - delttaP)) * Uc(+1)/Uc(+12);

    @#endif


@#endif 


@#elseif optimal == 2

kP = kPss;
AP = delttaP * kPss; 
GP = AP;
ap = AP/Y;
qP * (1 -  psi/2*((AP/AP(-1))-1)^2 -  psi*((AP/AP(-1))-1)*AP/AP(-1))
    +betta*Uc(+1)/Uc*qP(+1)* psi*((AP(+1)/AP)-1)*(AP(+1)/AP)^2 = 1;



@#elseif optimal == 3

        @#if time_to_build == 1
[name='Capital motion for protection capital']
AP = apss * Y;

[name='Public investment in protection capital']
GP = AP;

[name='Capital motion']
kP = (1 - delttaP) * kP(-1) + AP;


        @#elseif time_to_build == 4
[name='Capital motion for protection capital']
AP = apss * Y;

[name='Public investment in protection capital']
GP = AP(-3);

[name='Capital motion']
kP = (1 - delttaP) * kP(-1) + AP(-3);


        @#else
[name='Capital motion for protection capital']
AP = apss * Y;

[name='Public investment in protection capital']
GP = AP(-11);

[name='Capital motion']
kP = (1 - delttaP) * kP(-1) + AP(-11);

        @#endif

qP = 1;

ap = apss;



@#endif





%%%% ------- Natural disaster --------%%%%

[name='Damage function with protection capital']
damage = (1 + b1 * (kP(-1))^b2)/(1 + b1 * (kP(-1))^b2 + a1 * (temp - tempopt) + a2 * (temp - tempopt)^2);

[name='Climatic shock']
temp = (1 - rhot) * tempm + rhot * temp(-1)  + epsilon;


[name='Share of adaptative capital over total public capital']
share = kP/(kP + kI);


[name='Ratio of adaptation capital over GDP']
kp = kP/Y;


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

kI = kIss; 
AI = delttaI * kIss;
GI = AI;
qI = qIss;
ai = AI/Y;

kP = kPss;
AP = delttaP * kPss; 
GP = AP;
qP = qPss;
ap = AP/Y;

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

W = alfa * Y/n;
R = (1 - alfa) * Y/k;
c = C/Y;
in = Invest/Y;
share = kP/(kP + kI);
Rb = 1/betta - 1;
B = (mu * Yss - 2 * GI - 2 * GP)/(Rb - lambda);
b = B/Y;
T = lambda * B + mu * Yss - GI - GP;
s = (T - GI - GP)/Y;

kp = kP/Y;
t = T/Y;

S = s * Y;
G = GI + GP;
g = G/Y;
end; 








/* 
                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        V - SIMULATIONS
                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/


shocks;
var epsilon; stderr .01;
end;

//write_latex_dynamic_model;

resid(1);
steady; 
check;
options_.TeX=1;

stoch_simul(nomoments, nocorr, nofunctions, irf = 120) Y C n k Invest W R kI AI kP AP G B T S Rb  c in ai ap g b t s;
