function [ TMDCmatrix, TMDCmatrix_inv ] = TMDC_transformation_matrix( matrix_ind)
% Reference citation: Phys. Rev. B 92, 205108 (2015).
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017

% k space transformation: |k_QFT> = exp(i k r_i) |k_w>


sq2=sqrt(2);
sq3=sqrt(3);

%new basis:

noddx=4;
noddy=5;
noddz=3;
noddxz=1;
noddyz=2;
nevenx=10;
neveny=11;
nevenz=9;
nevenz2=6;
nevenxy=7;
nevenx2y2=8;

%old basis:
odz2=1;
odxy=2;
odx2y2=3;
odxz=4;
odyz=5;
os1x=6;
os1y=7;
os1z=8;
os2x=9;
os2y=10;
os2z=11;

mat_wan_oe=zeros(11,11);

mat_wan_oe(odxz,noddxz)=1.0;
mat_wan_oe(odyz,noddyz)=1.0;
mat_wan_oe(odz2,nevenz2)=1.0;
mat_wan_oe(odxy,nevenxy)=1.0;
mat_wan_oe(odx2y2,nevenx2y2)=1.0;

mat_wan_oe(os1x,noddx)=1/sq2;
mat_wan_oe(os2x,noddx)=-1/sq2;
mat_wan_oe(os1y,noddy)=1/sq2;
mat_wan_oe(os2y,noddy)=-1/sq2;
mat_wan_oe(os1z,noddz)=1/sq2;
mat_wan_oe(os2z,noddz)=1/sq2;

mat_wan_oe(os1x,nevenx)=1/sq2;
mat_wan_oe(os2x,nevenx)=1/sq2;
mat_wan_oe(os1y,neveny)=1/sq2;
mat_wan_oe(os2y,neveny)=1/sq2;
mat_wan_oe(os1z,nevenz)=1/sq2;
mat_wan_oe(os2z,nevenz)=-1/sq2;

% from px py pz basis to spherical harmonics basis
%L=1
ang_1=[-1/sq2,i/sq2,0;0,0,1;1/sq2,i/sq2,0];

%L=2
ang_2=[0,-i/sq2,1/sq2,0,0;0,0,0,-1/sq2,i/sq2;1,0,0,0,0;0,0,0,1/sq2,i/sq2;0,i/sq2,1/sq2,0,0];

%mirror symmetry xz plane in atomic basis
mirror_xz=diag([1,-1,1,1,-1,1,-1,1,1,-1,1]);
mirror_yz=diag([1,-1,1,-1,1,-1,1,1,-1,1,1]);

% in odd/even basis
rotate_pi=diag([-1,1,-1,-1,1,1,-1,1,1,1,-1]);

switch matrix_ind
    case 0
        % odd/even basis to wannier atomic basis
        TMDCmatrix = mat_wan_oe;
        TMDCmatrix_inv = inv(mat_wan_oe);
        
        
    case 1
        % odd/even basis to spherical basis
        TMDCmatrix = zeros(11);
        TMDCmatrix(1:5,1:5)=ang_2;
        TMDCmatrix(6:8,6:8)=ang_1;
        TMDCmatrix(9:11,9:11)=ang_1;
        TMDCmatrix=TMDCmatrix*mat_wan_oe;
        
        TMDCmatrix_inv=inv(TMDCmatrix);
        
        
    case 2
        % wannier basis to spherical basis
        TMDCmatrix = zeros(11);
        TMDCmatrix(1:5,1:5)=ang_2;
        TMDCmatrix(6:8,6:8)=ang_1;
        TMDCmatrix(9:11,9:11)=ang_1;
        
        TMDCmatrix_inv=inv(TMDCmatrix);
        
        
        
        
    case 3
        % rotation (120 degrees clockwise) in odd/even basis on states
        TMDCmatrix = zeros(11);
        TMDCmatrix(3,3)=1;
        TMDCmatrix(6,6)=1;
        TMDCmatrix(9,9)=1;
        TMDCmatrix(1,1)=-1/2;
        TMDCmatrix(2,1)=-sq3/2;
        TMDCmatrix(1,2)=sq3/2;
        TMDCmatrix(2,2)=-1/2;
        
        TMDCmatrix(4,4)=-1/2;
        TMDCmatrix(5,4)=-sq3/2;
        TMDCmatrix(4,5)=sq3/2;
        TMDCmatrix(5,5)=-1/2;
        
        TMDCmatrix(10,10)=-1/2;
        TMDCmatrix(11,10)=-sq3/2;
        TMDCmatrix(10,11)=sq3/2;
        TMDCmatrix(11,11)=-1/2;
        
        TMDCmatrix(7,7)=-1/2;
        TMDCmatrix(8,7)=-sq3/2;
        TMDCmatrix(7,8)=sq3/2;
        TMDCmatrix(8,8)=-1/2;
        
        TMDCmatrix_inv=inv(TMDCmatrix);
        
        
        
        
    case 4
        % rotation (60 degrees clockwise) in odd/even basis
        TMDCmatrix = zeros(11);
        TMDCmatrix(3,3)=1;
        TMDCmatrix(6,6)=1;
        TMDCmatrix(9,9)=1;
        
        TMDCmatrix(1,1)=1/2;
        TMDCmatrix(2,1)=-sq3/2;
        TMDCmatrix(1,2)=sq3/2;
        TMDCmatrix(2,2)=1/2;
        
        TMDCmatrix(4,4)=1/2;
        TMDCmatrix(5,4)=-sq3/2;
        TMDCmatrix(4,5)=sq3/2;
        TMDCmatrix(5,5)=1/2;
        
        TMDCmatrix(10,10)=1/2;
        TMDCmatrix(11,10)=-sq3/2;
        TMDCmatrix(10,11)=sq3/2;
        TMDCmatrix(11,11)=1/2;
        
        TMDCmatrix(7,7)=-1/2;
        TMDCmatrix(8,7)=sq3/2;
        TMDCmatrix(7,8)=-sq3/2;
        TMDCmatrix(8,8)=-1/2;
        
        TMDCmatrix_inv=inv(TMDCmatrix);
        
    case 5
        % rotation (60 degrees clockwise) in wannier basis
        TMDCmatrix = zeros(11);
        TMDCmatrix(1,1)=1;
        TMDCmatrix(8,8)=1;
        TMDCmatrix(11,11)=1;
        
        TMDCmatrix(4,4)=1/2;
        TMDCmatrix(5,4)=-sq3/2;
        TMDCmatrix(4,5)=sq3/2;
        TMDCmatrix(5,5)=1/2;
        
        TMDCmatrix(6,6)=1/2;
        TMDCmatrix(7,6)=-sq3/2;
        TMDCmatrix(6,7)=sq3/2;
        TMDCmatrix(7,7)=1/2;
        
        TMDCmatrix(9,9)=1/2;
        TMDCmatrix(10,9)=-sq3/2;
        TMDCmatrix(9,10)=sq3/2;
        TMDCmatrix(10,10)=1/2;
        
        TMDCmatrix(2,2)=-1/2;
        TMDCmatrix(3,2)=sq3/2;
        TMDCmatrix(2,3)=-sq3/2;
        TMDCmatrix(3,3)=-1/2;
        
        TMDCmatrix_inv=inv(TMDCmatrix);
        
    case 6
        %mirror symmetry x->x y->-y operation for atomic basis
        TMDCmatrix = mirror_xz;
        TMDCmatrix_inv=inv(TMDCmatrix);
        
    case 7
        %mirror symmetry x->x y->-y operation for odd/even basis
        TMDCmatrix = mat_wan_oe'*mirror_xz*mat_wan_oe;
        TMDCmatrix_inv=inv(TMDCmatrix);
        
    case 8
        %convert odd even X spin up down
        TMDCmatrix=zeros(22);
        TMDCmatrix(1:5,1:5)=eye(5);
        TMDCmatrix(17:22,6:11)=eye(6);
        TMDCmatrix(12:16,12:16)=eye(5);
        TMDCmatrix(6:11,17:22)=eye(6);
        
        TMDCmatrix_inv=inv(TMDCmatrix);
        
        
    case 9
        % rotation (120 degrees clockwise) in atomic basis on states
        TMDCmatrix = zeros(11);
        TMDCmatrix(3,3)=1;
        TMDCmatrix(6,6)=1;
        TMDCmatrix(9,9)=1;
        TMDCmatrix(1,1)=-1/2;
        TMDCmatrix(2,1)=-sq3/2;
        TMDCmatrix(1,2)=sq3/2;
        TMDCmatrix(2,2)=-1/2;
        
        TMDCmatrix(4,4)=-1/2;
        TMDCmatrix(5,4)=-sq3/2;
        TMDCmatrix(4,5)=sq3/2;
        TMDCmatrix(5,5)=-1/2;
        
        TMDCmatrix(10,10)=-1/2;
        TMDCmatrix(11,10)=-sq3/2;
        TMDCmatrix(10,11)=sq3/2;
        TMDCmatrix(11,11)=-1/2;
        
        TMDCmatrix(7,7)=-1/2;
        TMDCmatrix(8,7)=-sq3/2;
        TMDCmatrix(7,8)=sq3/2;
        TMDCmatrix(8,8)=-1/2;
        
        TMDCmatrix=mat_wan_oe*TMDCmatrix*inv(mat_wan_oe);
        
        TMDCmatrix_inv=inv(TMDCmatrix);
        
        
        
    case 10
        % mirror symmetry x->x y->-y operation for atomic basis
        TMDCmatrix = mirror_yz;
        TMDCmatrix_inv=inv(TMDCmatrix);
        
    case 11
        % mirror symmetry x->x y->-y operation for odd/even basis
        TMDCmatrix = mat_wan_oe'*mirror_yz*mat_wan_oe;
        TMDCmatrix_inv=inv(TMDCmatrix);
        
    case 12
        % rotate around x axis y->-y, z->-z for atomic basis
        TMDCmatrix= mat_wan_oe*rotate_pi*mat_wan_oe';
        TMDCmatrix_inv=inv(TMDCmatrix);
        
    case 13
        % rotate around x axis y->-y, z->-z for odd/even basis
        TMDCmatrix = rotate_pi;
        TMDCmatrix_inv=inv(TMDCmatrix);
        
        
    otherwise
        TMDCmatrix=0;
        TMDCmatrix_inv=0;
        
        
end


end

