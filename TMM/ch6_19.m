clearvars
clc

G = 0.8e11; % in N/m^2
kg=1000;
gr=2; % Gear ratio

nele = 3;   % no. of elements
d = [0.02 0 0.015];  % diameters of shafts in m

         
    connect = [1 1 2;   %% First Column is element number
               2 2 3
               3 3 4];    % Second & Third Column are Nodes (in sequence)
               
           
     coord = [1 0.00;       % first column node second column x coodinate
              2 0.40;
              3 0.40
              4 0.80];
         
     mass = [0 3  2 0];     % mass of disc at each node

     Dia = [0 0.2  0.1 0];  % dia of disc at node
     
 bc = 2;         %1: free-free ; 2:fixed-fixed ; 3:fixed-free ; 4:free-fixed
 for i=1:nele+1
    Ip(i) = mass(i)*Dia(i)^2/8;   % polar moment of inertia of each disc
end
 
 M = zeros(nele+1,nele+1);   %initializing global mass and stiffness matrix
 K = zeros(nele+1,nele+1);

 for i = 1:nele              %loop for finding elemental mass &
   
    nd1 = connect(i,2);     %stiffness matrix and assembly
    nd2 = connect(i,3);
    x1 = coord(nd1,2);
    x2 = coord(nd2,2);
    l = x2-x1;
    J = pi*d(i)^4/32;
    if(i==1)
        kele(:,:,i) = G*J/l*[1 -1;-1 1];
    elseif (i==2)
        kele(:,:,i) = kg*[1 -1;-1 1];  
    else
        kele(:,:,i) = G*(J/l)*[1 -1; -1 1];   % geared shaft
    end
    if(i==1)
        mele(:,:,i) = [Ip(nd1) 0; 0 Ip(nd2)];  % discs are included both sides
    elseif(i==2)% geared shaft
         mele(:,:,i) = [0 0; 0 0];
    else
        mele(:,:,i) = [Ip(nd1) 0;0 Ip(nd2)];  %disc is considered in both side
    end
    
    %assembly
    vec = [nd1 , nd2];  %global dof vector for assembly
    for ii = 1:2
        for jj = 1:2
            K(vec(ii),vec(jj)) =K(vec(ii),vec(jj))+kele(ii,jj,i);
            M(vec(ii),vec(jj)) =M(vec(ii),vec(jj))+mele(ii,jj,i);
        end
    end
 end
% imposing boundary condition

if(bc==1)       %free-free
    Kred = K;
    Mred = M;
elseif(bc==2)      %fixed-fixed
    Kred = K(2:nele,2:nele);
    Mred = M(2:nele,2:nele);
elseif(bc==3)       %fixed-free
    Kred = K(2:nele+1,2:nele+1);
    Mred = M(2:nele+1,2:nele+1);
else                %free-fixed
    Kred = K(1:nele,1:nele);
    Mred = M(1:nele,1:nele);
end
D = Mred\Kred;
[eig_vec,eig_val] = eig(D);

for i=1:2
    r=eig_val(i,i);
    if(r<0)
        wnf(i)=0;
    else
    wnf(i) = sqrt(r);   %critical speed
    end
end
%post processing

disp('Solution is printed to a text file "Exercise6.19_output_FEM.txt"');
fid = fopen('Output_Exercise6.19_FEM.txt','w');
fprintf(fid,'Finite Element method\n\n');    
fprintf(fid,'Number of elements = %d\n',nele);

if(nele <= 6)
for i=1:nele
    fprintf(fid,'Element-%d\n',i);
    fprintf(fid,'----------\n\n');
    fprintf(fid,'Mass Matrix [M]%d\n',i);
    for ii = 1:2
        for jj=1:2
            fprintf(fid,'%f\t',mele(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Stiffness matrix [K]%d\n',i);
    for ii = 1:2
        for jj=1:2
            fprintf(fid,'%f\t',kele(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end


fprintf(fid,'Global Mass matrix \n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf(fid,'%f \t',M(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n');
fprintf(fid,'Global Stiffness matrix \n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf(fid,'%f \t',K(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
if(bc==1)
    fprintf(fid, '"Free-Free" Boundary condition\n');
elseif(bc==2)      
   fprintf(fid, 'Fixed-Fixed Boundary condition\n');
elseif(bc==3)      
   fprintf(fid, 'Fixed-Free Boundary condition\n');
else                
    fprintf(fid, 'Free-Fixed Boundary condition\n');
end
fprintf(fid,'\n');
fprintf(fid,'Natural frequencies:\n');
fprintf(fid,'%.3f \n',wnf);
fprintf(fid,'\n');
fclose(fid);
end

