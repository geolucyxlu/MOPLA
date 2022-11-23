function [E_bar,M_bar,Em,Sm,e,s,vis,vis0,ind,idx] = MT_frictional_viscous_StressBC(Sigma,a,q,eta,eta0,Nc,Nm,n,phi,ys)
% MORI_TANAKA_FRICTION_VISCOUS_STRESSBC 
% 08/July/2021
%  description: To get bulk rehology for composite material with a connected
%  matrix using Mori-Tanaka homogenization appraoch. This model considers
%  the competition between frictional and viscous deformation 
%  Stress boundary condition
%  Input Parameters:-------------------------------------------------------
% Sigma:Imposed macroscale stress tensor, a 3-by-3 matrix.
% n:	The number of inclusions.
% a:	The initial three semi-axes of clasts (a1>a2>a3; a3=1), a 3-by-3 matrix.
% q:	The initial orientations of clasts as transformation matrix, 
%       a 3-by-3-by-n matrix. 
% eta:	Clasts' effective viscosities defined at the initial state, 
%       such as the imposed macroscopic strain-rate state, a 1-by-n vector.
% eta0: Matrix viscosity defined at the imposed satate.
% Nc:	Clasts' stress exponents, a 1-by-n vector.
% Nm:	Matrix stress exponent.
% phi:  Clast volume fraction.
% Output Parameters:-------------------------------------------------------
% E_bar: bulk strain rate
% M_bar: bulk compliance tensor
% Em:    Matrix strain rate
% Sm:    Matrix stress
% e:     strain rates in clasts
% s:     stresses in clasts
% vis:   clasts' viscosities at the convergent state
% vis0:  matrix viscosity at the convergent state
% -------------------------------------------------------------------------
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~] = FourIdentity();
   
   % matrix stiffness
   Cm = 2*eta0*Jd;
   Mm = FourTensorInv(Cm);
   %  strain rate invariants at which the viscosities of RDEs are defiend  
   REF= Inva(Sigma)*ones(1,n); 
   REF0 =Inva(Sigma);

   Ce = zeros(3,3,3,3,n);
   Me = Ce;
  for j=1:n
      Ce(:,:,:,:,j) = 2*eta(j)*Jd;
      Me(:,:,:,:,j) = 1/(2*eta(j))*Jd;
  end
  B_dil = Ce;
  mB    = Ce;
  s  = zeros(3,3,n); 
  e = s;
  ind=zeros(1,n);
  sI =zeros(1,n);
  idx=0;
  alpha=ind;
 %% Mori-Tanaka
 for ii=1:1000
     for kk=1:1000         
         for r=1:n
              SS=SnPI_vis(a(:,r));
              %  transform the Eshelby tensor(SS, PI) to the global coordinate
              S=Transform(SS,q(:,:,r)');
              t1 = Contract(S,Contract(Mm,Ce(:,:,:,:,r))); %secant compliances
              t2 = Jd-S+Nm*t1;
              t3 = FourTensorInv(t2);
              t4 = Jd+(Nm-1)*S;                                                  
              A_dil = Contract(t3,t4);     %A in Eq.23 Jiang,2014
              t5 = Contract(Ce(:,:,:,:,r),A_dil);
              B_dil(:,:,:,:,r)= Contract(t5,Mm);
         end
            BB = (1-phi)*Jd + phi*mean(B_dil,5);
            Bm_MT = FourTensorInv(BB);
            Sm = Multiply(Bm_MT,Sigma);
         for r=1:n
             Br_MT= Contract(B_dil(:,:,:,:,r) ,Bm_MT);
             s(:,:,r) = Multiply(Br_MT,Sigma);
             sI(r)= Inva(s(:,:,r));
             % switch rheology based on the local stress
             if ind(r)==1
                 ev = sI(r)/eta(r)/2;
                 zz = ys/ev/2;
             else
                 zz = (sI(r)/REF(r))^(1-Nc(r))*eta(r);
                 %REF(r)   = sI(r);
             end 
%              % check viscosity
%              alpha(r) = abs(zz/eta(r) - 1);
             eta(r)   = zz;
             %check stress
             alpha(r)    = abs(sI(r)/REF(r) - 1);
             REF(r)   = sI(r);
             Ce(:,:,:,:,r)  = 2*eta(r)*Jd;
         end

         SmI = Inva(Sm);
         % switch rheology based on the local stress
         if idx==1
             ev0= SmI/eta0/2;
             yy = ys/ev0/2;
         else
             yy = (SmI/REF0)^(1-Nm)*eta0;
             %REF0 = SmI;
         end
%          % check viscosity
%          beta = abs(yy/eta0 - 1);
         eta0    = yy;   
         % check stress
         beta    = abs(SmI/REF0 - 1);
         REF0    = SmI;
         Cm  = 2*eta0*Jd; 
         Mm = FourTensorInv(Cm);

         if ~any(alpha>0.001) && ~(beta>0.001)
             break
         end
         
         if kk==1000
              if any(alpha>0.001)&& beta < 0.001
                 disp('Warning: Clast not converged')
                 disp(find(alpha>0.001))
              elseif ~any(alpha>0.001) && beta > 0.001
                 disp('Warning: Matrix not converged')
                 YY=['beta = ',num2str(beta)];
                 disp(YY) 
              elseif any(alpha>0.001) && beta > 0.001
                 disp('Warning: Clast and Matrix not converged')
                 YY=['beta = ',num2str(beta)];
                 disp(YY) 
                 disp(find(alpha>0.001))
              end
          end 
  
     end
     
     for r=1:n
         Me(:,:,:,:,r)  = 1/(2*eta(r))*Jd;
         e(:,:,r) = Multiply(Me(:,:,:,:,r),s(:,:,r));
         mB(:,:,:,:,r)=Contract(Me(:,:,:,:,r),B_dil(:,:,:,:,r));
     end
     mB1 = (1-phi)*Mm + phi*mean(mB,5);
     M_bar = Contract(mB1,Bm_MT);
     E_bar = Multiply(M_bar,Sigma);
     Em = Multiply(Mm,Sm);
    
     % check whether local stress exceeds yield strength 
     idx = SmI>=ys;
     for k=1:n
        if sI(1,k)>=ys
            ind(1,k)=1;
        end
     end
      % make sure both clast and matrix stresses are caped by yielding strength  
     if ~any(sI>ys+0.01) && ~(SmI>ys+0.01)
        break
     end
     
end   
    vis=eta;
    vis0=eta0;
end

