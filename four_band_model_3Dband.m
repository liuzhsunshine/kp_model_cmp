%liuzhao@ustc
%04-05-2018

close all;
clear all;
clc;

%%%% dimension of the system %%%%%%%
N=4;

%%%%Pauli matrix
%%%%sub for sublattice degree of freedom%%%%%%
%%%%pes for pseudospin or orbital degree of freedom%%%%
%%%%val for valley degree of freedom%%%%%%%
%%%%spi for spin declgree of freedom%%%%%%%%%
unitm=eye(2,2);
sub_x=[0, 1;
       1, 0];
sub_y=[0, -i;
        i, 0];
sub_z=[1, 0;
       0,-1];
val_x=[0, 1;
       1, 0];
val_y=[0, -i;
        i, 0];
val_z=[1, 0;
       0,-1];

%%%% efficient of Pauli matrix %%%%%%
a13=1;
a20=1;
a30=0.2;
a33=0.1;
%%%% manifold range %%%%%%
stepx=0.01;
stepy=0.01;
kx=-0.5:stepx:0.5;
ky=-0.5:stepy:0.5;

d=zeros(N,length(kx),length(ky));

%%%% calculation the berry curvature%%%
for r=1:length(kx)
    for j=1:length(ky)
        %%%%Total Hamiltonian
        H=a13*kx(r)*kron(sub_x,val_z)+a20*ky(j)*kron(sub_y,unitm)+a30*kron(sub_z,unitm)+a33*kron(sub_z,val_z);
 
        %%%%Solve total Hamiltonian            
        d(:,r,j)=real(eig(H));
    end
end

d=sort(real(d));    %here v_val, d_sort, d_index are same size 

d1(:,:)=d(1,:,:);
d2(:,:)=d(2,:,:);
d3(:,:)=d(3,:,:);
d4(:,:)=d(4,:,:);
%%% plot 3D band structure
mesh(kx,ky,d1);
hold on;
mesh(kx,ky,d2);
hold on;
mesh(kx,ky,d3);
hold on;
mesh(kx,ky,d4);
grid on;
xlabel('kx');
ylabel('ky');
zlabel('E');
axis equal;


                    