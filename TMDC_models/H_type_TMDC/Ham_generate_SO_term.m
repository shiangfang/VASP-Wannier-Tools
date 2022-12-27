function [ Hmat ] = Ham_generate_SO_term( soc_m,soc_x )
% Reference citation: Phys. Rev. B 92, 205108 (2015).
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017


%   generate SOC term in odd/even basis

sq2=sqrt(2);
sq3=sqrt(3);
sq5=sqrt(5);
sq6=sqrt(6);
sq1=sqrt(10);

Hmat=zeros(22);

% rearrange states into M spin up; M spin down; X spin up... and so on,
Mmat2=zeros(22);
Mmat2(1:5,1:5)=eye(5);
Mmat2(6:10,12:16)=eye(5);
Mmat2(11:13,6:8)=eye(3);
Mmat2(14:16,17:19)=eye(3);
Mmat2(17:19,9:11)=eye(3);
Mmat2(20:22,20:22)=eye(3);
[Umat,Umatinv]= TMDC_transformation_matrix(0);

% from odd/even basis to atomic basis
Mmat3=zeros(22);
Mmat3(1:11,1:11)=Umat(:,:);
Mmat3(12:22,12:22)=Umat(:,:);

Mmat1=Mmat2*Mmat3;

MSO=diag([soc_m*ones(1,6),-1.5*soc_m*ones(1,4),0.5*soc_x*ones(1,4),-soc_x*ones(1,2),0.5*soc_x*ones(1,4),-soc_x*ones(1,2)]);

T1mat=zeros(10);
T2mat=zeros(6);

T1mat=[-1/sq2,i/sq2,0,0,0,0;
    0,0,sq2/sq3,-1/sq6,i/sq6,0;
    1/sq6,i/sq6,0,0,0,sq2/sq3;
    0,0,0,1/sq2,i/sq2,0;
    0,0,-1/sq3,-1/sq3,i/sq3,0;
    -1/sq3,-i/sq3,0,0,0,1/sq3;];


T2mat=[0,-i/sq2,1/sq2,0,0,0,0,0,0,0;
    0,0,0,-sq2/sq5,i*sq2/sq5,0,-i/sq1,1/sq1,0,0;
    sq3/sq5,0,0,0,0,0,0,0,-1/sq5,i/sq5;
    0,0,0,1/sq5,i/sq5,sq3/sq5,0,0,0,0;
    0,i/sq1,1/sq1,0,0,0,0,0,sq2/sq5,i*sq2/sq5;
    0,0,0,0,0,0,i/sq2,1/sq2,0,0;
    0,0,0,1/sq1,-i/sq1,0,-i*sq2/sq5,sq2/sq5,0,0;
    -sq2/sq5,0,0,0,0,0,0,0,-sq3/sq1,i*sq3/sq1;
    0,0,0,-sq3/sq1,-i*sq3/sq1,sq2/sq5,0,0,0,0;
    0,-i*sq2/sq5,-sq2/sq5,0,0,0,0,0,1/sq1,i/sq1];



Tmat=zeros(22);
Tmat(1:10,1:10)=T2mat;
Tmat(11:16,11:16)=T1mat;
Tmat(17:22,17:22)=T1mat;

Mattmp=Tmat*Mmat1;

Hmat=(Mattmp')*MSO*Mattmp;

%max(abs(Tmat'*Tmat-eye(22,22)))
%max(abs(Mmat1'*Mmat1-eye(22)))


end

