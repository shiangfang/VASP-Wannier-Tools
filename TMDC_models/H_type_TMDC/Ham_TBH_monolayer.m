function [ Hmat ] = Ham_TBH_monolayer( know,tbparms,a1,a2 )
% Reference citation: Phys. Rev. B 92, 205108 (2015).
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017

%take the symmetric minimu tb parms set
sq3=sqrt(3);
sq2=sqrt(2);

internal_dim=11;
Hmat=zeros(internal_dim);

% as vectors defined in the paper

del1=a1;
del2=(a1+a2);
del3=a2;

del4=-(2*a1+a2)/3;
del5=(a1+2*a2)/3;
del6=(a1-a2)/3;

del7=-2*(a1+2*a2)/3;
del8=2*(2*a1+a2)/3;
del9=2*(a2-a1)/3;

% all irreducible parameters

ep1=tbparms(1);
ep3=tbparms(2);
ep4=tbparms(3);
ep6=tbparms(4);
ep7=tbparms(5);
ep9=tbparms(6);
ep10=tbparms(7);
t111=tbparms(8);
t122=tbparms(9);
t133=tbparms(10);
t144=tbparms(11);
t155=tbparms(12);
t166=tbparms(13);
t177=tbparms(14);
t188=tbparms(15);
t199=tbparms(16);
t11010=tbparms(17);
t11111=tbparms(18);
t135=tbparms(19);
t168=tbparms(20);
t1911=tbparms(21);
t112=tbparms(22);
t134=tbparms(23);
t145=tbparms(24);
t167=tbparms(25);
t178=tbparms(26);
t1910=tbparms(27);
t11011=tbparms(28);
t541=tbparms(29);
t532=tbparms(30);
t552=tbparms(31);
t596=tbparms(32);
t5116=tbparms(33);
t5107=tbparms(34);
t598=tbparms(35);
t5118=tbparms(36);
t696=tbparms(37);
t6116=tbparms(38);
t698=tbparms(39);
t6118=tbparms(40);


% generate all parameters

ep2=ep1;
ep5=ep4;
ep8=ep7;
ep11=ep10;

% X-X, M-M type

t211=0.25*t111+0.75*t122;
t244=0.25*t144+0.75*t155;
t277=0.25*t177+0.75*t188;
t21010=0.25*t11010+0.75*t11111;

t222=0.75*t111+0.25*t122;
t255=0.75*t144+0.25*t155;
t288=0.75*t177+0.25*t188;
t21111=0.75*t11010+0.25*t11111;

t233=t133;
t266=t166;
t299=t199;

t235=(sq3/2)*t134-0.5*t135;
t268=(sq3/2)*t167-0.5*t168;
t2911=(sq3/2)*t1910-0.5*t1911;

t335=-(sq3/2)*t134-0.5*t135;
t368=-(sq3/2)*t167-0.5*t168;
t3911=-(sq3/2)*t1910-0.5*t1911;

t212=(sq3/4)*(t111-t122)-t112;
t245=(sq3/4)*(t144-t155)-t145;
t278=(sq3/4)*(t177-t188)-t178;
t21011=(sq3/4)*(t11010-t11111)-t11011;

t312=-(sq3/4)*(t111-t122)-t112;
t345=-(sq3/4)*(t144-t155)-t145;
t378=-(sq3/4)*(t177-t188)-t178;
t31011=-(sq3/4)*(t11010-t11111)-t11011;

t234=0.5*t134+(sq3/2)*t135;
t267=0.5*t167+(sq3/2)*t168;
t2910=0.5*t1910+(sq3/2)*t1911;

t334=0.5*t134-(sq3/2)*t135;
t367=0.5*t167-(sq3/2)*t168;
t3910=0.5*t1910-(sq3/2)*t1911;

% X-M type
t441=0.25*t541+0.75*t552;
t4107=0.25*t5107+0.75*t5118;

t452=0.75*t541+0.25*t552;
t4118=0.75*t5107+0.25*t5118;

t451=-(sq3/4)*t541+(sq3/4)*t552;
t442=t451;
t4117=-(sq3/4)*t5107+(sq3/4)*t5118;
t4108=t4117;

t431=-(sq3/2)*t532;
t497=-(sq3/2)*t598;

t432=-0.5*t532;
t498=-0.5*t598;

t496=t596;
t4106=-sq3*t5116/2;
t4116=-t5116/2;

Hmat=zeros(11);
Hmat_ext=zeros(11);

%%%%%%%%%%%%%%

% type X-X M-M

tbh_type1=[3,5,t135,t235,t335;
    6,8,t168,t268,t368;
    9,11,t1911,t2911,t3911];

for inds=1:3
    Hmat(tbh_type1(inds,1),tbh_type1(inds,2))=2*tbh_type1(inds,3)*cos(dot(know,del1))+tbh_type1(inds,4)*(exp(-i*dot(know,del2))+exp(-i*dot(know,del3)))+tbh_type1(inds,5)*(exp(i*dot(know,del2))+exp(i*dot(know,del3)));
end

tbh_type2=[1,2,t112,t212,t312;
    3,4,t134,t234,t334;
    4,5,t145,t245,t345;
    6,7,t167,t267,t367;
    7,8,t178,t278,t378;
    9,10,t1910,t2910,t3910;
    10,11,t11011,t21011,t31011];

for inds=1:7
    Hmat(tbh_type2(inds,1),tbh_type2(inds,2))=-2*i*tbh_type2(inds,3)*sin(dot(know,del1))+tbh_type2(inds,4)*(exp(-i*dot(know,del2))-exp(-i*dot(know,del3)))+tbh_type2(inds,5)*(-exp(i*dot(know,del2))+exp(i*dot(know,del3)));
end


% type X-M
tbh_type3=[3,1,t431;
    5,1,t451;
    4,2,t442;
    10,6,t4106;
    9,7,t497;
    11,7,t4117;
    10,8,t4108];

for inds=1:7
    Hmat(tbh_type3(inds,1),tbh_type3(inds,2))=tbh_type3(inds,3)*(exp(i*dot(know,del4))-exp(i*dot(know,del6)));
end
    
tbh_type4=[4,1,t441,t541;
    3,2,t432,t532;
    5,2,t452,t552;
    9,6,t496,t596;
    11,6,t4116,t5116;
    10,7,t4107,t5107;
    9,8,t498,t598;
    11,8,t4118,t5118];

for inds=1:8
    Hmat(tbh_type4(inds,1),tbh_type4(inds,2))=tbh_type4(inds,3)*(exp(i*dot(know,del4))+exp(i*dot(know,del6)))+tbh_type4(inds,4)*exp(i*dot(know,del5));
end


%%%%%%%%%%%%%%
Hmat_ext(9,6)=t696*(exp(i*dot(know,del7))+exp(i*dot(know,del8))+exp(i*dot(know,del9)));
Hmat_ext(11,6)=t6116*(exp(i*dot(know,del7))-0.5*exp(i*dot(know,del8))-0.5*exp(i*dot(know,del9)));
Hmat_ext(10,6)=(0.5*sq3)*t6116*(-exp(i*dot(know,del8))+exp(i*dot(know,del9)));

Hmat_ext(9,8)=t698*(exp(i*dot(know,del7))-exp(i*dot(know,del8))/2-exp(i*dot(know,del9))/2);
Hmat_ext(9,7)=0.5*sq3*t698*(0-exp(i*dot(know,del8))+exp(i*dot(know,del9)));
Hmat_ext(10,7)=(3/4)*t6118*(0+exp(i*dot(know,del8))+exp(i*dot(know,del9)));
Hmat_ext(11,7)=(sq3/4)*t6118*(0+exp(i*dot(know,del8))-exp(i*dot(know,del9)));
Hmat_ext(10,8)=(sq3/4)*t6118*(0+exp(i*dot(know,del8))-exp(i*dot(know,del9)));
Hmat_ext(11,8)=t6118*(exp(i*dot(know,del7))+exp(i*dot(know,del8))/4+exp(i*dot(know,del9))/4);

Hmat=Hmat+Hmat'+Hmat_ext+Hmat_ext';


tbh_diag=[1,1,ep1,t111,t211;
    2,2,ep2,t122,t222;
    3,3,ep3,t133,t233;
    4,4,ep4,t144,t244;
    5,5,ep5,t155,t255;
    6,6,ep6,t166,t266;
    7,7,ep7,t177,t277;
    8,8,ep8,t188,t288;
    9,9,ep9,t199,t299;
    10,10,ep10,t11010,t21010;
    11,11,ep11,t11111,t21111];

for inds=1:11
   Hmat(tbh_diag(inds,1),tbh_diag(inds,2))=tbh_diag(inds,3)+2*tbh_diag(inds,4)*cos(dot(know,del1))+2*tbh_diag(inds,5)*(cos(dot(know,del2))+cos(dot(know,del3)));
end



end

