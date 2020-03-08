classdef powergrid_state_estimation
   % Power Grid State Estimation
   %
   % The essence of the problem of estimate powergrid state is to make an 
   % unbiased estimation of the variables of equations when the number of 
   % equations is greater than the number of variables.
   %
   % For powergrid systems, the variable is the node voltage (i.e.,the
   % state value of the powergrid, which consists of the real part and the 
   % imaginary part), and part of the node and branch power (containing 
   % errors due to objective conditions) are the quantity measurement, that
   % will be given. The functional structure of the system of equations can 
   % be derived from Kirchhoff's theorem.
   % 
   % This code provides two solutions to solve the above problem, namely 
   % the least square method and the fast decoupled method. In the least
   % square method, the Newton iteration is adopt to iteratively calculate 
   % the node voltage. In the meantime, the fast decoupled method
   % simplified and decoupled the functions of least square method
   % according to the actual characteristic of high-voltage powergrid
   % system.
   %
   % The input of this code is the modified standard IEEE 30-bus data, and
   % the filename is 'iSE30Bus.txt'. We manually adjusted the input values 
   % of node and branch power to simulate the errors contained in the data 
   % under real conditions. After running over the code, the estimated 
   % powergrid data will be output to 'oStateEstimation.txt'.
   %
   % Copyright (c) 2020 Yaoyu Chen
   
   
   properties
       num_Node; num_Branc;  % Number of nodes (bus) and branches
       num_kBus; num_kBranc; % Number of given node power and brance power
       PQ; Branc;            % Given node power and brance power
       G; B;                 % Real and imaginary matrix of node admittance
       Res; V;               % Branch impedance and node voltage
       F; J;                 % inferred power values and Jacobian matrix  
   end
   
   
   methods
       
       % Shared functions
       % Initial function
       function self = powergrid_state_estimation(path)
           self = ReadPGStat(self, path);
           self = CalcAdmitMat(self);
           self.V = [ones(self.num_Node,1),zeros(self.num_Node,1)];
       end
       
       
       % Infering the node voltage
       function infer(self, alg, iter_num, thrs, save_path)
           % --------------------------------------------------------------
           % Two methods are given to solve the problem of estimate 
           % powergrid states, one is the classical least-square method and 
           % the other is fast-decoupled methd. After calculate the bus 
           % voltage iteratively, the inferred data (e.g., node voltage, 
           % node power, branch power and branch current, etc) will be 
           % saved according to the given filename.
       
           if isempty(iter_num)
               iter_num = 100;
           end
           if isempty(thrs)
               thrs = 1e-3;
           end
       
           if strcmp(alg, 'least-square')
               self = least_square(self, iter_num, thrs);
           elseif strcmp(alg, 'fast-decoupled')
               self = fast_decoupled(self, iter_num, thrs);
           else
               error('The chioce of algorithm is not expected \n');
           end
           SavePGStat(self, save_path);
       end
       
       
       % Reading powergrid data
       function self = ReadPGStat(self, path)
           % --------------------------------------------------------------
           % The standard IEEE 30-Bus system data are ultilized in this 
           % example.
           % Parameters in the file:
           %    num_Node, num_Branc: number of nodes, branches
           %    num_kBus:            number of given node power
           %    num_kBranch:         number of given branch power
           %    Res:                 impedances of branches
           %    PQ, Branc:           node power and branch power
       
           if isempty(path)
              error('The filename of powergrid data must be given! \n')
           end
               
           fprintf('reading data...\n'); fp = fopen(path, 'r');
           fscanf(fp,'%s',1);
           self.num_Node   = fscanf(fp,'%d',1);
           self.num_Branc  = fscanf(fp,'%d',1);
           self.num_kBus   = fscanf(fp,'%d',1);
           self.num_kBranc = fscanf(fp,'%d',1);

           fscanf(fp,'%s',1);  % read branch parameters (impedances)
           self.Res = zeros(self.num_Branc,4);
           for i = 1:self.num_Branc
               for j = 1:4
                 self.Res(i,j) = fscanf(fp,'%f',1);
               end
           end

           fscanf(fp,'%s',1);  % read the given bus parameters (node power)
           self.PQ = zeros(self.num_kBus,3);
           for i = 1:self.num_kBus
               for j = 1:3
                   self.PQ(i,j) = fscanf(fp,'%f',1);
               end
           end
           sortrows(self.PQ);

           fscanf(fp,'%s',1);  % read the given branch power
           self.Branc = zeros(self.num_kBranc,4);
           for i = 1:self.num_kBranc
               for j = 1:4
                 self.Branc(i,j) = fscanf(fp,'%f',1);
               end
           end
           fclose(fp);
       end
       
       
       % Calculate the node admittance matrix
       function self = CalcAdmitMat(self)
           % --------------------------------------------------------------
           % Calculate the node admittance matrix, which will be used both 
           % in the state estimation equations and its derivative equations.
           % In this function, G and B represent the real and imaginary 
           % part of the node admittance matrix respectively.
          
           G = zeros(self.num_Node,self.num_Node);
           B = zeros(self.num_Node,self.num_Node);  %#ok<*PROPLC>
           Res = self.Res;
           for k = 1:self.num_Branc
               i =  Res(k,1); j = Res(k,2);
               f =  Res(k,3)/(Res(k,3)^2 + Res(k,4)^2); % branch conductance
               g = -Res(k,4)/(Res(k,3)^2 + Res(k,4)^2); % branch susceptance

               G(i,i) = G(i,i) + f; G(j,j) = G(j,j) + f;
               B(i,i) = B(i,i) + g; B(j,j) = B(j,j) + g;
               G(i,j) = G(i,j) - f; G(j,i) = G(j,i) - f;
               B(i,j) = B(i,j) - g; B(j,i) = B(j,i) - g;
           end
           self.G = G; self.B = B;
       end
       
       
       % Save the calculated results of the power grid state
       function SavePGStat(self, save_path)
           if isempty(save_path)
               save_path = 'results.txt';
           end
           
           fp = fopen(save_path,'w'); V = self.V;
           fprintf(fp,'\n\n        *****  State Estimation Results  *****\n\n');
           fprintf(fp,'     Bus        Vamp     Ang(Deg)    P           Q\n');
           for i = 1:self.num_Node
               f = sqrt(V(i,1)^2 + V(i,2)^2);
               g = atan(V(i,2)/V(i,1))*180/pi;
               fprintf(fp,'%8d%12.6f%13.6f%12.6f%12.6f\n',i,f,g,self.F(2*i-1),self.F(2*i));
           end

           % Output branch power flow
           mwLt = 0.0; mvLt = 0.0;
           fprintf(fp,'\n\n       i       j      Pij         Qij         Pji...         Qji       MW Loss     MVar Loss    Iij        Iji\n');
           for k = 1:self.num_Branc
               i = self.Res(k,1); j = self.Res(k,2);
               r = self.Res(k,3); x = self.Res(k,4);

               f   = complex(V(i,1),V(i,2));  % Vi
               gi  = (complex(V(i,1),V(i,2)) - complex(V(j,1),V(j,2)))/complex(r,x);  % Iij
               pij = real(f*conj(gi)); qij = imag(f*conj(gi));
               f   = complex(V(j,1),V(j,2));  % Vj
               gj  = (complex(V(j,1),V(j,2)) - complex(V(i,1),V(i,2)))/complex(r,x);  % Iji
               pji = real(f*conj(gj)); qji = imag(f*conj(gj));

               mwL  = pij + pji;  mvL  = qij + qji;
               mwLt = mwLt + mwL; mvLt = mvLt+mvL;
               fprintf(fp,'%8d%8d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n',i,j,pij,qij,pji,qji,mwL,mvL,gi,gj);
           end

           fprintf(fp,'\n     Total Loss:       %12.6f  +j%12.6f',mwLt,mvLt);
           fclose(fp);
       end
       
       
       % Functions of least-square method
       % Inferred the powergrid state using least-square method
       function self = least_square(self, iter_num, thrs)
           % --------------------------------------------------------------
           % In the least-square method, the Newton iteration is adopt
           % to iteratively calculate the node voltage V. During the 
           % iteration, the node and branch power can be estimated base on 
           % the current V, then the difference between the estimated power
           % value and the ground truth can be further ultilized to modify V.

           iter = 0;
           while iter < iter_num
               self = CalcStatEsEq(self);
               self = CalcJacb(self);

               gt_F = [self.PQ(:,2:3);self.Branc(:,3:4)]'; % ground truth value
               MSErr = pinv(self.J'*self.J) * self.J' * (gt_F(:) - self.F);
               err = max(abs(MSErr));
               if err < thrs
                   break;
               else
                   self.V = self.V + reshape(MSErr,2,self.num_Node)';
                   iter = iter + 1;
               end
           end
           if iter >= iter_num
               error('Error! Iteration Divergence! \n');
           end
       end
       
       
       % Calculate node and branch power using state estimation equations
       function self = CalcStatEsEq(self)
           % --------------------------------------------------------------
           % Calculate the inferred node and branch power using node 
           % admittance matrix and iterative node voltage.
           % By comparing the calculated values with given truth values, 
           % the result of approximation of the iterative node voltage can 
           % be evaluated.
           
           self.F = zeros(2*(self.num_kBus + self.num_kBranc),1);
           G = self.G; B = self.B; V = self.V; Res = self.Res; count = 0;
           for n = 1:self.num_kBus   % Node power
               i = self.PQ(n,1);
               f = G(i,:)*V(:,1) - B(i,:)*V(:,2);
               g = G(i,:)*V(:,2) + B(i,:)*V(:,1);

               self.F(count+1) = V(i,1)*f + V(i,2)*g;
               self.F(count+2) = V(i,2)*f - V(i,1)*g;
               count = count + 2;
           end

           for k = 1:self.num_kBranc  % Branch power
               i = self.Branc(k,1); j = self.Branc(k,2); flag = 0;
               for n = 1:self.num_Branc
                   n1 = Res(n,1); n2 = Res(n,2);
                   if (i==n1 && j==n2) || (i==n2 && j==n1)
                       r = Res(n,3)^2 + Res(n,4)^2;
                       f = (V(i,1)-V(j,1))*Res(n,3) + (V(i,2)-V(j,2))*Res(n,4);
                       g = (V(i,1)-V(j,1))*Res(n,4) - (V(i,2)-V(j,2))*Res(n,3);

                       self.F(count+1) = (V(i,1)*f - V(i,2)*g)/r;
                       self.F(count+2) = (V(i,2)*f + V(i,1)*g)/r;
                       count = count + 2; flag = 1; break;
                   end
               end
               if flag == 0
                   error('Error in calculate F! \n');
               end
           end
       end
       
       
       % Calculate Jacobian matrix
       function self = CalcJacb(self)
           % --------------------------------------------------------------
           % Calculate Jacobian matrix of the state estimation equations.
           % The Jacobian matrix are ultilized in the gradient calculation
           % of Newton iteration. In this function, the number of rows of 
           % the Jacobian matrix corresponds to the number of state 
           % estimation equations, while the number of columns corresponds 
           % to the number of unknown variables (ie, node voltage).
       
           self.J = zeros(2*(self.num_kBus + self.num_kBranc), 2*self.num_Node);
           i = 0; Res = self.Res; G = self.G; B = self.B; V = self.V;
           
           for k = 1:self.num_kBus
               j = 0; m = self.PQ(k,1);
               for n = 1:self.num_Node
                   if m == n
                       g = G(m,:)*V(:,1) - B(m,:)*V(:,2);
                       f = G(m,:)*V(:,2) + B(m,:)*V(:,1);
                   else
                       g = 0; f = 0;
                   end
                   self.J(i+1,j+1) =  G(m,n)*V(m,1) + B(m,n)*V(m,2) + g;  % dP/de
                   self.J(i+1,j+2) =  G(m,n)*V(m,2) - B(m,n)*V(m,1) + f;  % dP/df
                   self.J(i+2,j+1) =  G(m,n)*V(m,2) - B(m,n)*V(m,1) - f;  % dQ/de
                   self.J(i+2,j+2) = -G(m,n)*V(m,1) - B(m,n)*V(m,2) + g;  % dQ/df
                   j = j + 2;
                end
                i = i + 2;
           end

           for m = 1:self.num_kBranc
               i = self.Branc(m,1); j = self.Branc(m,2); flag = 0;
               for n = 1:self.num_Branc
                   n3 = Res(n,1); n4 = Res(n,2);

                   if (i==n3 && j==n4) || (j==n3 && i==n4)
                       pos = 2*(self.num_kBus + m); flag = 1;
                       r = Res(n,3)^2 + Res(n,4)^2;
                       f = V(i,1)-V(j,1); g = V(i,2)-V(j,2);

                       self.J(pos-1,i*2-1) = ( V(i,1)*Res(n,3) - V(i,2)*Res(n,4) +...
                           Res(n,3)*f + Res(n,4)*g)/r;  %dPij/dei
                       self.J(pos-1,j*2-1) = (-V(i,1)*Res(n,3) + V(i,2)*Res(n,4))/r; %dPij/dej
                       self.J(pos-1,i*2)   = ( V(i,1)*Res(n,4) + V(i,2)*Res(n,3) +...
                           Res(n,3)*g - Res(n,4)*f)/r;  %dPij/dfi
                       self.J(pos-1,j*2)   = (-V(i,1)*Res(n,4) - V(i,2)*Res(n,3))/r; %dPij/dfj

                       self.J(pos,i*2-1) = ( V(i,1)*Res(n,4) + V(i,2)*Res(n,3) -...
                           Res(n,3)*g + Res(n,4)*f)/r;  %dQij/dei
                       self.J(pos,j*2-1) = (-V(i,1)*Res(n,4) - V(i,2)*Res(n,3))/r;  %dQij/dej
                       self.J(pos,i*2)   = (-V(i,1)*Res(n,3) + V(i,2)*Res(n,4) +...
                           Res(n,3)*f + Res(n,4)*g)/r;  %dQij/dfi
                       self.J(pos,j*2)   = ( V(i,1)*Res(n,3) - V(i,2)*Res(n,4))/r;  %dQij/dfj
                       break;
                   end
               end
               if flag == 0
                   error('Error in calculate J(branch)! \n');
               end
           end
       end
       
       
       % Functions of fast-decoupled method
       % Inferred the powergrid state using fast-decoupled method
       function self = fast_decoupled(self, iter_num, thrs)
           % --------------------------------------------------------------
           % In the actual high-voltage power system, it can be found that 
           % the real and imaginary parts of the node voltage have quite 
           % different influence on the real and imaginary parts of the 
           % power. Meanwhile, the phase angle difference of the voltage 
           % across the node is quite small. Hence, with some simple 
           % assumptions, the complexity of the powergrid state equation 
           % can be greatly reduced, and the state equations (i.e.,node 
           % and branch power) can also be decoupled (into real and 
           % imaginary parts) and calculated separately.
           
           [Ba, Br] = CalcBaBr(self); iter = 0;
           while iter < iter_num
               self = CalcFDEsEq(self);
               [delV, err] = CalcErr(self, Ba, Br);

               if err(1)<=thrs && err(2)<=thrs
                   break
               end
               if err(1) > thrs
                   self.V(:,2) = self.V(:,2) + delV(:,1);
               end
               if err(2) > thrs
                   self.V(:,1) = self.V(:,1) + delV(:,2);
               end
               iter = iter + 1;
           end
           if iter >= iter_num
               error('Error! Iteration Divergence! \n');
           end
       end
       
       
       % Calculate the decoupled Jacobian matrix
       function [Ba, Br] = CalcBaBr(self)
           % --------------------------------------------------------------
           % In Fast-decoupled method, the Jacobian matrix can also be 
           % decoupled into two part, and they are both constant matrices
           % .i.e., J = | Ba  0 |
           %            | 0  Br | 
           
           Res = self.Res;
           Ba = zeros(self.num_kBus+self.num_kBranc, self.num_Node);
           Ba(1:self.num_kBus, 1:self.num_Node) = -self.B(self.PQ(:,1),:);
           Br = zeros(self.num_kBus+self.num_kBranc, self.num_Node);
           Br(1:self.num_kBus, 1:self.num_Node) = -self.B(self.PQ(:,1),:);

           for k = 1:self.num_kBranc
               i = self.Branc(k,1); j = self.Branc(k,2); flag = 0;
               for n = 1:self.num_Branc
                   n1 = Res(n,1); n2 = Res(n,2);
                   if (i==n1 && j==n2) || (i==n2 && j==n1)
                       r = -Res(n,4)/(Res(n,3)^2 + Res(n,4)^2);
                       Ba(self.num_kBus+k,i) = -r;
                       Ba(self.num_kBus+k,j) = r;
                       flag = 1;
                   end
               end
               if flag == 0
                  error('Error in calculat Ba! \n');
               end
           end

           for k = 1:self.num_kBranc
               i = self.Branc(k,1); j = self.Branc(k,2); flag = 0;
               for n = 1:self.num_Branc
                   n1 = Res(n,1); n2 = Res(n,2);
                   if (i==n1 && j==n2) || (i==n2 && j==n1)
                       r = -Res(n,4)/(Res(n,3)^2 + Res(n,4)^2);
                       Br(self.num_kBus+k,i) = -r;
                       Br(self.num_kBus+k,j) = r;
                       flag = 1;
                   end
               end
               if flag == 0
                  error('Error in calculat Br! \n');
               end
           end
       end
       
       
       % Calculate the node and branch power in fast-decoupled method
       function self = CalcFDEsEq(self)
           % --------------------------------------------------------------
           % In the fast-decoupled method, the state estimation equations
           % are simplified a lot. Thus, the equation to estimate node
           % and branch power need to adjust into the following form.
           
           self.F = zeros(2*(self.num_kBus + self.num_kBranc), 1);
           V = self.V; G = self.G; B = self.B; Res = self.Res;
           for count = 1:self.num_kBus
               i = self.PQ(count,1); f = 0; g = 0;
               for j = 1:self.num_Node
                   m = sin(V(i,2)-V(j,2)); n = cos(V(i,2)-V(j,2));
                   f = f + V(i,1) * V(j,1) * (G(i,j)*n + B(i,j)*m);
                   g = g + V(i,1) * V(j,1) * (G(i,j)*m - B(i,j)*n);
               end
               self.F(count) = f;
               self.F(self.num_kBus + self.num_kBranc + count) = g;
           end

           for count = 1:self.num_kBranc
               i = self.Branc(count,1); j = self.Branc(count,2); flag = 0;
               for n = 1:self.num_Branc
                   n1 = Res(n,1); n2 = Res(n,2);

                   if i == n1 && j == n2
                       f = Res(n,3)/(Res(n,3)^2 + Res(n,4)^2);
                       g = -Res(n,4)/(Res(n,3)^2 + Res(n,4)^2); flag = 1;

                       self.F(self.num_kBus+count) = V(i,1)*V(i,1)*f -...
                           V(i,1)*V(j,1)*(f*cos(V(i,2)-V(j,2)) +...
                           g*sin(V(i,2)-V(j,2)));
                       self.F(2*self.num_kBus+self.num_kBranc+count) =...
                           -V(i,1)*V(i,1)*g -V(i,1)*V(j,1)*...
                           (f*sin(V(i,2)-V(j,2))-g*cos(V(i,2)-V(j,2)));         

                   elseif i==n2 && j==n1
                       f = Res(n,3)/(Res(n,3)^2 + Res(n,4)^2);
                       g = -Res(n,4)/(Res(n,3)^2 + Res(n,4)^2); flag = 1;

                       self.F(self.num_kBus+count) = V(i,1)*V(i,1)*f +...
                           V(i,1)*V(j,1)*(-f*cos(V(j,2)-V(i,2)) +...
                           g*sin(V(j,2)-V(i,2)));
                       self.F(2*self.num_kBus+self.num_kBranc+count) =...
                           -V(i,1)*V(i,1)*g +V(i,1)*V(j,1)*...
                           (f*sin(V(j,2)-V(i,2))+g*cos(V(j,2)-V(i,2)));
                   end
               end
               if flag == 0
                  error('Error in calculate F!');
               end
           end
       end
       
       
       % Calculate the error of the node voltage
       function [delV, err] = CalcErr(self, Ba, Br)
           % --------------------------------------------------------------
           % Compare the difference estimated power and the ground truth
           % value, and then calculated the values that need to be 
           % corrected for the node voltage.
           
           delV = zeros(self.num_Node,2); err = zeros(1,2);
           Err1 = [self.PQ(:,2);self.Branc(:,3)] -...
               self.F(1:self.num_kBus+self.num_kBranc);
           Err2 = [self.PQ(:,3);self.Branc(1:self.num_kBranc,4)] -...
               self.F(self.num_kBus+self.num_kBranc+1:end);

           delV(:,1) = pinv(Ba'*Ba)*Ba'*Err1; err(1) = max(abs(delV(:,1)));
           delV(:,2) = pinv(Br'*Br)*Br'*Err2; err(2) = max(abs(delV(:,2)));
       end
   end
end