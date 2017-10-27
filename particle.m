% Upgraded from tinp.m to study disorder and shell
clear all;
close all;
clc;

% Set parameters and initialise matrices
L=4; % should be odd for a sphere, s.t. there is a lattice point in the centre
task='sphere'; % default cube, specify 'sphere' otherwise
R=2;  % only used if task='sphere' 

% Conversion factor between angstrom and lattice points
fac=36; m0=-0.28; m2=32/(fac*fac); A=3/fac; Rcutsq=1.5; % to find nearest neighbors
N=L*L*L;
% initialising
coords=zeros(N,3);
incoords=zeros(N,3); 
H=zeros(4*N,4*N);

% Set lattice points on square grid, e.g. for L=2 the order is
% (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1)
for i=1:L
    for j=1:L
       for k=1:L
            incoords(L*L*(i-1)+L*(j-1)+k,1)=i-1;
            incoords(L*L*(i-1)+L*(j-1)+k,2)=j-1;
            incoords(L*L*(i-1)+L*(j-1)+k,3)=k-1;
            coords(L*L*(i-1)+L*(j-1)+k,1)=i-1;
            coords(L*L*(i-1)+L*(j-1)+k,2)=j-1;
            coords(L*L*(i-1)+L*(j-1)+k,3)=k-1;
       end
    end 
end

% Count only atoms within the sphere 
if task=='sphere'
    N_sphere=0;
    for i=1:N
        if ((incoords(i,1)-(L-1)/2)^2+(incoords(i,2)-(L-1)/2)^2+(incoords(i,3)-(L-1)/2)^2)<R*R
            N_sphere=N_sphere+1;
        end
    end
    coords=zeros(N_sphere,3);
    dummy=1; % dummy counts atoms within
    for i=1:N
        if ((incoords(i,1)-(L-1)/2)^2+(incoords(i,2)-(L-1)/2)^2+(incoords(i,3)-(L-1)/2)^2)<R*R
            coords(dummy,1)=incoords(i,1);
            coords(dummy,2)=incoords(i,2);
            coords(dummy,3)=incoords(i,3);
            dummy=dummy+1;
        end
    end
    N=N_sphere;
    H=zeros(4*N,4*N);
    incoords=zeros(N,3);
    dummy=1;
    for i=1:N
        incoords(dummy,1)=coords(i,1);
        incoords(dummy,2)=coords(i,2);
        incoords(dummy,3)=coords(i,3);
        dummy=dummy+1;
    end
end

% Calculate Hamiltonian
for i=1:N
    for j=1:N
        if j==i
            H(4*(j-1)+1,4*(j-1)+1)=(m0+6*m2);
            H(4*(j-1)+2,4*(j-1)+2)=-(m0+6*m2);
            H(4*(j-1)+3,4*(j-1)+3)=(m0+6*m2);
            H(4*(j-1)+4,4*(j-1)+4)=-(m0+6*m2);
        else
            if ((coords(i,1)-coords(j,1)).^2+(coords(i,2)-coords(j,2)).^2+(coords(i,3)-coords(j,3)).^2) < Rcutsq
                H(4*(i-1)+1,4*(j-1)+1)=-m2;
                H(4*(i-1)+2,4*(j-1)+2)=m2;
                H(4*(i-1)+3,4*(j-1)+3)=-m2;
                H(4*(i-1)+4,4*(j-1)+4)=m2;
                                
                if (coords(i,1)>coords(j,1) || coords(i,2)>coords(j,2) || coords(i,3)>coords(j,3))
                    H(4*(i-1)+1,4*(j-1)+4)=-1i*A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=-1i*A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=-1i*A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=-1i*A/2;
                    
                    H(4*(i-1)+1,4*(j-1)+4)=H(4*(i-1)+1,4*(j-1)+4)-A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=H(4*(i-1)+2,4*(j-1)+3)-A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=H(4*(i-1)+3,4*(j-1)+2)+A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=H(4*(i-1)+4,4*(j-1)+1)+A/2;
                    
                    H(4*(i-1)+1,4*(j-1)+2)=-1i*A/2;
                    H(4*(i-1)+2,4*(j-1)+1)=-1i*A/2;
                    H(4*(i-1)+3,4*(j-1)+4)=1i*A/2;
                    H(4*(i-1)+4,4*(j-1)+3)=1i*A/2;        
                elseif (coords(i,1)<coords(j,1) || coords(i,2)<coords(j,2) || coords(i,3)<coords(j,3))
                    H(4*(i-1)+1,4*(j-1)+4)=1i*A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=1i*A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=1i*A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=1i*A/2;
                  
                    H(4*(i-1)+1,4*(j-1)+4)=H(4*(i-1)+1,4*(j-1)+4)+A/2;
                    H(4*(i-1)+2,4*(j-1)+3)=H(4*(i-1)+2,4*(j-1)+3)+A/2;
                    H(4*(i-1)+3,4*(j-1)+2)=H(4*(i-1)+3,4*(j-1)+2)-A/2;
                    H(4*(i-1)+4,4*(j-1)+1)=H(4*(i-1)+4,4*(j-1)+1)-A/2;
                    
                    H(4*(i-1)+1,4*(j-1)+2)=1i*A/2;
                    H(4*(i-1)+2,4*(j-1)+1)=1i*A/2;
                    H(4*(i-1)+3,4*(j-1)+4)=-1i*A/2;
                    H(4*(i-1)+4,4*(j-1)+3)=-1i*A/2;
                else
                end
            end
        end
    end
end

% diagonalise Hamiltonian and pick eigenvalues close to zero
[eigvecs,eigvals]=eig(H);
    
for i=1:4*N
    if eigvals(i,i)<0 && eigvals(i+1,i+1)>0
        dummy=i+1;
    end
end
surf_states_tot=50; % how many states above/below Dirac point to pick
surf_eigvals=zeros(2*surf_states_tot,1);

for i=-surf_states_tot:surf_states_tot-1
        surf_eigvals(surf_states_tot+1+i,1)=eigvals(dummy+i,dummy+i);
end

% calculate radius from volume; output surface gap + analytical result
Reff=nthroot(3*N/(4*pi),3); 
fprintf('%d %d %f %f %f \n', N, Reff, eigvals(dummy,dummy), A/Reff, eigvals(dummy+2,dummy+2)-eigvals(dummy,dummy));

% calculate surface charge density due to a particular surface state
state_num=dummy+1;
surf_density=zeros(N,1);

for i=1:N
    surf_density(i,1)=eigvecs(4*(i-1)+1,state_num).*conj(eigvecs(4*(i-1)+1,state_num));
    surf_density(i,1)=surf_density(i,1)+eigvecs(4*(i-1)+2,state_num).*conj(eigvecs(4*(i-1)+2,state_num));
    surf_density(i,1)=surf_density(i,1)+eigvecs(4*(i-1)+3,state_num).*conj(eigvecs(4*(i-1)+3,state_num));
    surf_density(i,1)=surf_density(i,1)+eigvecs(4*(i-1)+4,state_num).*conj(eigvecs(4*(i-1)+4,state_num));
end

% plot eigenavalues and probability density for one eigenvector
set(gcf, 'Visible', 'off')
myfig=figure;
subplot(1,2,1)
plot(surf_eigvals,'o') 
hold on
plot(surf_states_tot-2*N+state_num,surf_eigvals(surf_states_tot-2*N+state_num),'ro')
hold off
plottitle_str=sprintf('State # %d/%d', state_num, 4*N);
title(plottitle_str)
ylabel('E (eV)');
ax = gca;
ax.XTickLabel = {};
subplot(1,2,2)
scatter3(coords(:,1)*fac/10,coords(:,2)*fac/10,coords(:,3)*fac/10, [ ], surf_density(:),'.')
colorbar('southoutside')
colormap(jet)
caxis([0,8./N])
xlabel('(nm)');
ylabel('(nm)');
zlabel('(nm)');
axis([-1*fac/10,(L+1)*fac/10,-1*fac/10,(L+1)*fac/10,-1*fac/10,(L+1)*fac/10]);
axis square;    
