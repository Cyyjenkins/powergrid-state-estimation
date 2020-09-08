%% Power Grid State Estimation
clc;    % Clear windows
format
%% Read power flow data from document "iSE30Bus.txt"
fp = fopen('iSE30Bus.txt','r');
fscanf(fp,'%s',1);          % read one sting
node = fscanf(fp,'%d',1);   % read number of nodes
nzl = fscanf(fp,'%d',1);    % read number of branches
kBus = fscanf(fp,'%d',1);   % read number of known Bus

fscanf(fp,'%s',1);          % read a string
zl = zeros(nzl,4);          % read branch parameters
for i = 1:nzl
    for j = 1:4
      zl(i,j) = fscanf(fp,'%f',1);
    end                     
end

fscanf(fp,'%s',1);          % read one string
N = zeros(node,3);          % read node powers
for i = 1:node   
    for j = 1:3
        N(i,j) = fscanf(fp,'%f',1);   % read generator parameters
    end
end
sortrows(N);
PQ = zeros(2*node,1);       % reshape the node power
for i = 1:node
    PQ(2*i-1) = N(i,2);
    PQ(2*i) = N(i,3);
end

fscanf(fp,'%s',1);          % read one string
Bus = zeros(kBus,4);        % save known Bus: i  j  Pij  Qij
for i = 1:kBus
    for j = 1:4
      Bus(i,j) = fscanf(fp,'%f',1);
    end                     % read generator parameters
end
fclose(fp);

%% Calculate the node admittance matrix
G = zeros(node,node);
B = zeros(node,node);
for k = 1:nzl
    i = zl(k,1);
    j = zl(k,2);
    h = zl(k,3)*zl(k,3) + zl(k,4)*zl(k,4);
    f = zl(k,3)/h;               % branch conductance
    g = -zl(k,4)/h;              % branch susceptance
    G(i,i) = G(i,i)+f;
    G(j,j) = G(j,j)+f;
    B(i,i) = B(i,i)+g;
    B(j,j) = B(j,j)+g;
    G(i,j) = G(i,j)-f;
    G(j,i) = G(j,i)-f;
    B(i,j) = B(i,j)-g;
    B(j,i) = B(j,i)-g;
end
G;                                % show
B;                                % show

%% Assign initial values for bus voltage e1 f1 e2 f2...
V = zeros(2*node,1);              % save: e1 f1 e2 f2 ...
for i = 1:node
   V(2*i-1) = 1;
end
V(59) = 0.9344;
V(60) = -0.2133;

%% State Estimation iteration
iter = 0;               % iteration number
F = zeros(2*node-2+2*kBus,1); % Generate the Power Flow Matrix z=h(x)
while iter<100                      
    k = 0;
    for i = 1:node-1    % Node Power            
            f = 0.0; g = 0.0;
            for j = 1:node
                f = f + G(i,j)*V(2*j-1) - B(i,j)*V(2*j); 
                g = g + G(i,j)*V(2*j) + B(i,j)*V(2*j-1); 
            end
            F(k+1) = V(2*i-1)*f + V(2*i)*g;
            F(k+2) = V(2*i)*f - V(2*i-1)*g;
            k = k+2;           
    end   
       
    for i = 1:kBus       % Branch Power
        j = Bus(i,1);m = Bus(i,2);
        for o = 1:nzl
            o1 = zl(o,1);o2 = zl(o,2);
            if j == o1
                if m == o2
                z1 = zl(o,3)*zl(o,3) + zl(o,4)*zl(o,4);
                d1 = V(2*j-1) - V(2*m-1);      % e(i)-e(j)
                d2 = V(2*j) - V(2*m);          % f(i)-f(j)
                d3 = d1*zl(o,3) + d2*zl(o,4);
                d4 = d1*zl(o,4) - d2*zl(o,3);
                d5 = V(2*j-1)*d3 - V(2*j)*d4;
                d6 = V(2*j)*d3 + V(2*j-1)*d4;
                
                F(k+1) = d5/z1;
                F(k+2) = d6/z1;
                k = k+2;
                end
            end
                
            if j == o2 
                if m == o1
                z1 = zl(o,3)*zl(o,3) + zl(o,4)*zl(o,4);
                d1 = V(2*j-1) - V(2*m-1);      % e(i)-e(j)
                d2 = V(2*j) - V(2*m);          % f(i)-f(j)
                d3 = d1*zl(o,3) + d2*zl(o,4);
                d4 = d1*zl(o,4) - d2*zl(o,3);
                d5 = V(2*j-1)*d3 - V(2*j)*d4;
                d6 = V(2*j)*d3 + V(2*j-1)*d4;
                
                F(k+1) = d5/z1;
                F(k+2) = d6/z1;
                k = k+2;   
                end
            end
        end     
    end 
   
       % Calculate Jacobian matrix: J=H(x) 
    J = zeros(2*node+2*kBus-2,2*node-2);
    i1 = 0;
    k = 0;
    for i = 1:node-1
            j1 = 0;
            for j = 1:node-1
                   if i ~= j
                       J(i1+1,j1+1) = G(i,j)*V(2*i-1)+B(i,j)*V(2*i);  % P-e
                       J(i1+1,j1+2) = G(i,j)*V(2*i)-B(i,j)*V(2*i-1);  % P-f
                       J(i1+2,j1+1) = G(i,j)*V(2*i)-B(i,j)*V(2*i-1);  % Q-e
                       J(i1+2,j1+2) = -G(i,j)*V(2*i-1)-B(i,j)*V(2*i); % Q-f
                   else
                       f = 0.0; g = 0.0;
                       for k = 1:node
                           f = f + G(i,k)*V(2*k-1) - B(i,k)*V(2*k);
                           g = g + G(i,k)*V(2*k) + B(i,k)*V(2*k-1);
                       end
                       J(i1+1,j1+1)=G(i,i)*V(2*i-1)+B(i,i)*V(2*i)+f;  % P-e
                       J(i1+1,j1+2)=G(i,i)*V(2*i)-B(i,i)*V(2*i-1)+g;  % P-f
                       J(i1+2,j1+1)=G(i,i)*V(2*i)-B(i,i)*V(2*i-1)-g;  % Q-e
                       J(i1+2,j1+2)=-G(i,i)*V(2*i-1)-B(i,i)*V(2*i)+f; % Q-f
                   end
                   j1 = j1+2;
             end
             i1 = i1+2;
    end

    for i = 1:kBus
        i6 = Bus(i,1);i2 = Bus(i,2);
        for i3 = 1:nzl     %find the impedance of the branch which p/q is given
            i4 = zl(i3,1);
            i5 = zl(i3,2);
            if i4 == i6
                if i5 == i2
                z = zl(i3,3)*zl(i3,3) + zl(i3,4)*zl(i3,4);
                i7 = V(2*i4-1) - V(2*i5-1);
                i8 = V(2*i4) - V(2*i5);
                
                J1=i7*zl(i3,3)+i8*zl(i3,4)+V(2*i4-1)*zl(i3,3)-V(2*i4)*zl(i3,4);
                J2=-V(2*i4-1)*zl(i3,3)+V(2*i4)*zl(i3,4);
                J3=V(2*i4-1)*zl(i3,4)+zl(i3,3)*i8-zl(i3,4)*i7+V(2*i4)*zl(i3,3);
                J4=-V(2*i4-1)*zl(i3,4)-V(2*i4)*zl(i3,3);
                J5=V(2*i4)*zl(i3,3)-zl(i3,3)*i8+zl(i3,4)*i7+V(2*i4-1)*zl(i3,4);
                J6=-V(2*i4)*zl(i3,3)-V(2*i4-1)*zl(i3,4);
                J7=i7*zl(i3,3)+i8*zl(i3,4)+V(2*i4)*zl(i3,4)-V(2*i4-1)*zl(i3,3);
                J8=-V(2*i4)*zl(i3,4)+V(2*i4-1)*zl(i3,3);
                               
                J(2*node-2+2*i-1,i4*2-1) = J1/z;  %Pij/ei
                J(2*node-2+2*i-1,i5*2-1) = J2/z;  %Pij/ej
                J(2*node-2+2*i-1,i4*2) = J3/z;    %Pij/fi
                J(2*node-2+2*i-1,i5*2) = J4/z;    %Pij/fj
                
                J(2*node-2+2*i,i4*2-1) = J5/z;    %Qij/ei
                J(2*node-2+2*i,i5*2-1) = J6/z;    %Qij/ej
                J(2*node-2+2*i,i4*2) = J7/z;      %Qij/fi
                J(2*node-2+2*i,i5*2) = J8/z;      %Qij/fj
                end
            end
        end
    end
  
    F1 = zeros(2*node-2+2*kBus,1);    %F1=PQ-F^(n)
    for i = 1:2*node-2
        F1(i) = PQ(i)-F(i);
    end
    for i =  1:kBus
        F1(2*node-2+2*i-1) = Bus(i,3)-F(2*node-2+2*i-1);
        F1(2*node-2+2*i) = Bus(i,4)-F(2*node-2+2*i);
    end
    
    X = zeros(2*node-2,1);            %derta X
    X1 = zeros(2*node-2,2*node-2);
    
    X1 = inv(J'*J);
    X = X1*J'*F1;
       
    Er = max(abs(X));
    if Er < 0.0000001
        break;     %  break the iterations if convergence
    else
        for i = 1:2*node-2
            V(i) = V(i) + X(i);
        end
        iter = iter+1;
    end
end

%% Output power flow results
fp=fopen('oStateEstimation.txt','w');
fprintf(fp,'\n\n                        *****  State Estimation Results  *****\n\n');
if iter>=100
    fprintf(fp,'    Iteration Divergence !    ');

else               
    % Output bus quantities
    fprintf(fp,'     Bus        Vamp     Ang(Deg)    P           Q\n');
    for i=1:node
        f=sqrt(V(2*i-1)*V(2*i-1)+V(2*i)*V(2*i));
        g=atan(V(2*i)/V(2*i-1))*180/3.1415926;
        fprintf(fp,'%8d%12.6f%13.6f%12.6f%12.6f\n',i,f,g,F(2*i-1),F(2*i));
    end

    % Output branch power flow
    mwLt=0.0; mvLt=0.0;
    fprintf(fp,'\n\n       i       j      Pij         Qij         Pji         Qji       MW Loss     MVar Loss    Iij        Iji\n');
    for k=1:nzl
        i=zl(k,1);
        j=zl(k,2);
        r=zl(k,3);
        x=zl(k,4);
        f=complex(V(2*i-1),V(2*i));      % Vi
        g1=(complex(V(2*i-1),V(2*i))-complex(V(2*j-1),V(2*j)))/complex(r,x);  % Iij
        pij=real(f*conj(g1));
        qij=imag(f*conj(g1));
        f=complex(V(2*j-1),V(2*j));      % Vj
        g2=(complex(V(2*j-1),V(2*j))-complex(V(2*i-1),V(2*i)))/complex(r,x);  % Iji
        pji=real(f*conj(g2));
        qji=imag(f*conj(g2));
        mwL=pij+pji;
        mvL=qij+qji;
        mwLt=mwLt+mwL;
        mvLt=mvLt+mvL;
        fprintf(fp,'%8d%8d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n',i,j,pij,qij,pji,qji,mwL,mvL,g1,g2);
    end
    
    fprintf(fp,'\n     Total Loss:       %12.6f  +j%12.6f',mwLt,mvLt);
    fprintf(fp,'\n     Iteration times:  %5d (times)\n',iter);
end
fclose(fp);
