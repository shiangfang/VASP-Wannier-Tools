% This code generates the tight-binding Hamiltonian and the corresponding
% band plots for four transition metal dichalcogenides (2H) in monolayer
% form without spin-orbit coupling
% 
% Reference citation: Phys. Rev. B 92, 205108 (2015).
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017

clear all;

% tmdc type variable(1:MoS2, 2:MoSe2, 3:WS2, 4:WSe2)
tmdc=1;
TMDC_monolayer;

%% band structure calculations

% Gamma-M-K-Gamma
k1=[0,0,0];
k2=b1/2;
k3=(2*b1-b2)/3;
knum=50;
scan_klist=[k1;k2;k3;k1];

[ all_kpts, scale_axis] = generate_k_line( knum, scan_klist );
knum_tot=size(all_kpts);
knum_tot=knum_tot(1);


tb_bands=zeros(11,knum_tot);
tb_vecs=zeros(11,11,knum_tot);
for ind = 1:knum_tot
    ind
    k_now=all_kpts(ind,:);
    
    [ Hmat_TB ] = Ham_TBH_monolayer( k_now,mono_tbh_parm,a1,a2 );
    
    
    [V,D]=eig(Hmat_TB);
    [tb_bands(:,ind),indband]=sort(real(diag(D)),'ascend');
    V=V(:,indband);
    tb_vecs(:,:,ind)=V(:,:);
    
end

%%

figure(1);
line([scale_axis(50),scale_axis(50)],[-100,100],'Color','k');
hold on
line([scale_axis(100),scale_axis(100)],[-100,100],'Color','k');

plot(scale_axis,tb_bands,'r','LineWidth',2);

%ylabel('Band Energy (eV)');
axis([-inf,inf,-7,5]);
box on;

ax = gca;
set(gca,'XTick',[0,scale_axis(50),scale_axis(100),1]);
set(gca,'YTick',[-6,-4,-2,0,2,4,6]);
set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
set(gca,'YTickLabel',{'-6','-4','-2','0','2','4','6'});
%set(gca,'YTickLabel',{'','','','','','',''});
set(gca,'FontSize',26);
set(gca,'Fontname','Times New Roman');
ylabel('E (eV)');


hold off;

