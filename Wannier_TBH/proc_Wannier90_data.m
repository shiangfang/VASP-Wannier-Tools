function [ r2mn,xyz_center,ham_r,ham_real,ham_imag,num_wann ] = proc_Wannier90_data()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

fname_xyz='wannier90_centres.xyz';
fname_r2mn='wannier90.r2mn';
fname_dat='wannier90_hr.dat';
fname_save='wannier90.DFT_KIT.mat';

% read r2mn data
data_r2mn_raw=dlmread(fname_r2mn);
data_r2mn=data_r2mn_raw(:,3);
num_wann=sqrt(length(data_r2mn));
r2mn=reshape(data_r2mn,num_wann,num_wann);

% process xyz data
fid1=fopen(fname_xyz,'r');
tmp1=fgets(fid1);
tmp1=fgets(fid1);
xyz_center=zeros(num_wann,3);
for inds=1:num_wann
    tmp1=fgets(fid1);
    tmp11=textscan(tmp1,'%s %f %f %f');
    xyz_center(inds,1)=tmp11{2};
    xyz_center(inds,2)=tmp11{3};
    xyz_center(inds,3)=tmp11{4};
end
fclose(fid1);


fid2=fopen(fname_dat,'r');
tmp1=fgets(fid2);
num_wann0=fscanf(fid2,'%d',[1]);
if num_wann0 ~= num_wann
    'ERROR'
end
num_Rs=fscanf(fid2,'%d',[1]);
R_degeneracy=fscanf(fid2,'%d',[num_Rs]);
% size(R_degeneracy)
ham_tmp=fscanf(fid2,'%f',[7,num_Rs*num_wann*num_wann]);
all_index=reshape(1:(num_Rs*num_wann*num_wann),num_wann*num_wann,num_Rs);

ham_r=zeros(num_Rs,3);
ham_imag=zeros(num_wann,num_wann,num_Rs);
ham_real=zeros(num_wann,num_wann,num_Rs);

for indR=1:num_Rs
    deg_fac=R_degeneracy(indR);
    tmpH=ham_tmp(:,all_index(:,indR));
    ham_r(indR,1:3)=[int32(tmpH(1)),int32(tmpH(2)),int32(tmpH(3))];
    tmp_real=tmpH(6,:)/deg_fac;
    tmp_imag=tmpH(7,:)/deg_fac;
    
    ham_real(:,:,indR)=reshape(tmp_real,num_wann,num_wann);
    ham_imag(:,:,indR)=reshape(tmp_imag,num_wann,num_wann);
    
end

fclose(fid2);

save(fname_save,'num_wann','r2mn','xyz_center','ham_r','ham_real','ham_imag');

end

