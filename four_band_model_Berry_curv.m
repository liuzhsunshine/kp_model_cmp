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
ds=stepx*stepy;      % area of squares

%%%% calculation the berry curvature%%%
omegaxy=zeros(N,length(kx),length(ky));
%omega1=zeros(length(kx),length(ky));
chern_sum=zeros(N,length(kx));
chern=zeros(N);
for r=1:length(kx)
    for j=1:length(ky)
        %%%%Total Hamiltonian
        H=a13*kx(r)*kron(sub_x,val_z)+a20*ky(j)*kron(sub_y,unitm)+a30*kron(sub_z,unitm)+a33*kron(sub_z,val_z);
        %fprintf('value of H is: ')  %test line, fprintf means print string
        %display(H)                  %test line, display here means print matrix
        %%%%Hx=partialH/partialkx
        %%%%Hy=partialH/partialky
        Hx=a13*kron(sub_x,val_z); 
        %fprintf('value of Hx is: ')  %test line
        %display(Hx)                  %test line
        Hy=a20*kron(sub_y,unitm);
        %fprintf('value of Hy is: ')  %test line
        %display(Hy)                  %test line
 
        %%%%Solve total Hamiltonian
        [v,d]=eig(H);               
        d_val=real(eig(H));
        [d_sort,d_index]=sort(d_val);    %here v_val, d_sort, d_index are same size
        v_sort=v(:,d_index);             %sort eigenvector(column) by order of eigenvalue
        %fprintf('value of d_val is: ')  %test line
        %display(d_val)                  %test line
        
        %%%%Normalize eigenvectors
        for t=1:N
            v_so(:,t)=v_sort(:,t)/norm(v_sort(:,t));
        end
        %fprintf('value of v_so is: ')   %test line
        %display(v_so)                   %test line
  
        %%%%Berry curvature of lower band  
        for m=1:N
            omega=0;   
            for n=1:N
                if n==m
                    omega=0;
                else
                %fprintf('value of n is: %d\n', n)    %test line
                %fprintf('value of m is: %d\n', m)    %test line      
                de=power(d_sort(m)-d_sort(n)+eps,-2); 
                %fprintf('value of d_sort(m) is: %d\n', d_sort(m))  %test line
                %fprintf('value of d_sort(n) is: %d\n', d_sort(n))  %test line
                %fprintf('value of de is: %d\n', de)                %test line
                mHxn=v_so(:,m)'*Hx*v_so(:,n);
                %fprintf('value of mHxn is: %d\n')        %test line
                %display(mHxn)                        %test line
                nHym=v_so(:,n)'*Hy*v_so(:,m);
                %fprintf('value of nHym is: %d\n')        %test line
                %display(nHym)                        %test line
                mHyn=v_so(:,m)'*Hy*v_so(:,n);
                %fprintf('value of mHyn is: %d\n')        %test line
                %display(mHyn)                        %test line
                nHxm=v_so(:,n)'*Hx*v_so(:,m);
                %fprintf('value of nHxm is: %d\n')        %test line
                %display(nHxm)                        %test line
                dH=mHxn*nHym-mHyn*nHxm;
                %fprintf('value of dH is: %d\n')      %test line
                %display(dH)                          %test line
                omega=omega+de*dH;                    %the value is pure imaginary
                %fprintf('value of omega is: %d\n')   %test line
                %display(omega)                       %test line
                end
            end
            omega(isnan(omega))=0;        %change NaN to zero
            omegaxy(m,r,j)=i*omega;   
            %fprintf('Berry curvature is: ')      %test line
            %display(omegaxy(m,r,j))              %test line
        end        
        chern_sum(1,r)=chern_sum(1,r)+power(3*pi,-1)*real(omegaxy(1,r,j))*ds;   
    end
    chern(1)=chern(1)+chern_sum(1,r);
end

%%% print Chern number of lower band
fprintf('Chern number of lower band is: %d\n', chern(1))

%%% plot the result
omega1(:,:)=real(omegaxy(1,:,:));       %Change complex to real to plot surface
mesh(kx,ky,omega1);
%surfl(kx,ky,omega1);
grid on;
xlabel('kx');
ylabel('ky');
zlabel('Berry curvature');

                    