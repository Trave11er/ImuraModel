clear all;
close all;
clc;

% Set parameters and initialise matrices
N=20;
kx=0.2; ky=0.0; % |kx| and |ky| should be less that Pi
% Conversion factor angstrom -> lattice points and other constants
m0=-1; m2=0.5; m1=0.5; A=2; B=2;
% initialising
H=zeros(4*N,4*N);
eigenvalues=zeros(4*N,1);

% Parameters for the bandstructure plot
alleigenvalues=zeros(4*N,10);
k_points=50;
k_spacing=0.0628;

for l=0:k_points-1
    kx=k_spacing*l;
    
    % Calculate Hamiltonian
    for i=1:N
        for j=1:N
            if j==i
                H(4*(i-1)+1,4*(j-1)+1)=(m0+2*m2*(3-cos(kx)-cos(ky)));
                H(4*(i-1)+2,4*(j-1)+2)=-(m0+2*m2*(3-cos(kx)-cos(ky)));
                H(4*(i-1)+3,4*(j-1)+3)=(m0+2*m2*(3-cos(kx)-cos(ky)));
                H(4*(i-1)+4,4*(j-1)+4)=-(m0+2*m2*(3-cos(kx)-cos(ky)));
                
                H(4*(i-1)+1,4*(j-1)+4)=A*(sin(kx)-1i*sin(ky));
                H(4*(i-1)+2,4*(j-1)+3)=A*(sin(kx)-1i*sin(ky));
                H(4*(i-1)+3,4*(j-1)+2)=A*(sin(kx)+1i*sin(ky));
                H(4*(i-1)+4,4*(j-1)+1)=A*(sin(kx)+1i*sin(ky));
            elseif j==i+1
                H(4*(i-1)+1,4*(j-1)+1)=-m1;
                H(4*(i-1)+2,4*(j-1)+2)=m1;
                H(4*(i-1)+3,4*(j-1)+3)=-m1;
                H(4*(i-1)+4,4*(j-1)+4)=m1;
                
                H(4*(i-1)+1,4*(j-1)+2)=-1i*B/2;
                H(4*(i-1)+2,4*(j-1)+1)=-1i*B/2;
                H(4*(i-1)+3,4*(j-1)+4)=1i*B/2;
                H(4*(i-1)+4,4*(j-1)+3)=1i*B/2;
            elseif j==i-1
                H(4*(i-1)+1,4*(j-1)+1)=-m1;
                H(4*(i-1)+2,4*(j-1)+2)=m1;
                H(4*(i-1)+3,4*(j-1)+3)=-m1;
                H(4*(i-1)+4,4*(j-1)+4)=m1;
                
                H(4*(i-1)+1,4*(j-1)+2)=1i*B/2;
                H(4*(i-1)+2,4*(j-1)+1)=1i*B/2;
                H(4*(i-1)+3,4*(j-1)+4)=-1i*B/2;
                H(4*(i-1)+4,4*(j-1)+3)=-1i*B/2;
            else
            end
        end
    end
    
    % diagonalise Hamiltonian and pick eigenvalues close to zero
    [eigvecs,eigvals]=eig(H);
    
    for i=1:4*N
        eigenvalues(i,1)=real(eigvals(i,i));
    end
    eigenvalues=sort(eigenvalues);
    alleigenvalues(:,l+1)=eigenvalues(:);
end

for i=0:k_points-1
    k=k_spacing*i;
    plot(k, alleigenvalues(:,i+1),'b.');
    hold on
    if k~=0
        plot(-k, alleigenvalues(:,i+1),'b.');
        hold on
    end
end
hold off

