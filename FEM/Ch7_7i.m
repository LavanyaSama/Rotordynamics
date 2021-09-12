clearvars
clc

rho=7800; % density 7800 kg/m^3
G = 0.8e11; % in N/m^2
len=3;
d = 0.03;   % diameters of shafts in m   
  n = [5,10,50,100];   % store 100 vectors
for loop=1:4
  nele = n(loop);  
connect = zeros(nele, 3);  % reserve memory for n [1x3] vectors
 coord = zeros(nele+1,2);
 coord(1,:)=[1 0];
for k = 1:nele
    v = [k k k+1];   % obtain the vector, your code comes here
    connect(k, :) = v;  % collect the vectors in a matrix
    coord(k+1,:)=[k+1 (k)*(len/nele)];
end
bc = 3;         %1: free-free ; 2:fixed-fixed ; 3:fixed-free ; 4:free-fixed

M = zeros(nele+1,nele+1);   %initializing global mass and stiffness matrix
K = zeros(nele+1,nele+1);

for i = 1:nele              %loop for finding elemental mass &
    nd1 = connect(i,2);     %stiffness matrix and assembly
    nd2 = connect(i,3);
    x1 = coord(nd1,2);
    x2 = coord(nd2,2);
     l = x2-x1;
    J = (pi*d^4)/32;
    kele(:,:,i) = G*J/l*[1 -1;-1 1];
    mele(:,:,i) =(rho*J*l/6)*[2 1; 1 2];
     
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

for i=1:nele
    r=eig_val(i,i);
    wnf(i) = sqrt(r); %critical speed
end
for p=1:length(wnf)-1
  for i =1:length(wnf)-1  %arranging in asceending order
    if wnf(i)>wnf(i+1)
        temp_wnf = wnf(i);
        wnf(i) = wnf(i+1);
        wnf(i+1) = temp_wnf;
        temp_mode(:,1) = mode(:,i);
        mode(:,i) = mode(:,i+1);
        mode(:,i+1) = temp_mode(:,1);
    end
  end
end
%post processing
fid = fopen('Output_Exercise7_7_FEM.txt','a');
fprintf(fid,'Finite Element method\n\n'); 
fprintf(fid,'\n');
fprintf(fid,'Fixed-Free Boundary Condition\n\n'); 
fprintf(fid,'Natural frequencies:\n');
fprintf(fid,'number of elements:\t');
fprintf(fid,'%d',nele);
fprintf(fid,'\n');
fprintf(fid,'%.3f \n',wnf);
fprintf(fid,'\n');
fclose(fid);
end