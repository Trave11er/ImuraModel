clear all;
close all;
clc;

% Set parameters and initialise matrices
L=12; M=16; % lengths along x and y (transverse) directions
k=0.0; % |k| should be less that Pi
m0=-1; B=0.7; m2=-m0; m1=m2;  A=1;
% initialising
N=L*M; 
Rcutsq=1.5; % to find nearest neighbors
H=zeros(4*N,4*N);
eigenvalues=zeros(4*N,1);
coords=zeros(N,2);

% Parameters for the bandstructure plot 
alleigenvalues=zeros(4*N,10);
k_points=40;
k_spacing=0.0157;
for l=0:k_points-1
    k=k_spacing*l;

% Set lattice points on square grid
for i=1:L
    for j=1:M
        coords(M*(i-1)+j,1)=i-1;
        coords(M*(i-1)+j,2)=j-1;
    end 
end

% Calculate Hamiltonian
for i=1:N
    for j=1:N
        if j==i
            H(4*(j-1)+1,4*(j-1)+1)=(m0+2*m1*(3-cos(k)));
            H(4*(j-1)+2,4*(j-1)+2)=-(m0+2*m1*(3-cos(k)));
            H(4*(j-1)+3,4*(j-1)+3)=(m0+2*m1*(3-cos(k)));
            H(4*(j-1)+4,4*(j-1)+4)=-(m0+2*m1*(3-cos(k)));
            
            H(4*(i-1)+1,4*(j-1)+2)=B*sin(k);
            H(4*(i-1)+2,4*(j-1)+1)=B*sin(k);
            H(4*(i-1)+3,4*(j-1)+4)=-B*sin(k);
            H(4*(i-1)+4,4*(j-1)+3)=-B*sin(k);
        else
            if ((coords(i,1)-coords(j,1)).^2+(coords(i,2)-coords(j,2)).^2) < Rcutsq
                H(4*(i-1)+1,4*(j-1)+1)=-m2;
                H(4*(i-1)+2,4*(j-1)+2)=m2;
                H(4*(i-1)+3,4*(j-1)+3)=-m2;
                H(4*(i-1)+4,4*(j-1)+4)=m2;
                                
                if (coords(i,1)>coords(j,1) || coords(i,2)>coords(j,2))
                    H(4*(i-1)+1,4*(j-1)+4)=-1i*A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=-1i*A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=-1i*A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=-1i*A/2;
                    
                    H(4*(i-1)+1,4*(j-1)+4)=H(4*(i-1)+1,4*(j-1)+4)-A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=H(4*(i-1)+2,4*(j-1)+3)-A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=H(4*(i-1)+3,4*(j-1)+2)+A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=H(4*(i-1)+4,4*(j-1)+1)+A/2;     
                elseif (coords(i,1)<coords(j,1) || coords(i,2)<coords(j,2))
                    H(4*(i-1)+1,4*(j-1)+4)=1i*A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=1i*A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=1i*A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=1i*A/2;
                  
                    H(4*(i-1)+1,4*(j-1)+4)=H(4*(i-1)+1,4*(j-1)+4)+A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=H(4*(i-1)+2,4*(j-1)+3)+A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=H(4*(i-1)+3,4*(j-1)+2)-A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=H(4*(i-1)+4,4*(j-1)+1)-A/2;
                else
                end
            end
        end
    end
end

% diagonalise Hamiltonian and pick eigenvalues close to zero
[eigvecs,eigvals]=eig(H);

for i=1:4*N
    eigenvalues(i,1)=real(eigvals(i,i));
    if imag(eigvals(i,i)>1e-8)
        display('Warning')
    end
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
%scatter(x, eigenvalues(:),'o');
ylim([-1,1]);
xlim([-0.628,0.628]);