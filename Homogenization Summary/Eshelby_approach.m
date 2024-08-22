function[S_bar,C_bar,e,s,vis]=Eshelby_approach(E,n,a,q,eta,eta0,Nk,Nm,phi)
%  Eshelby_approach.m       
%  description: To get bulk rehology and local fields using Eshelby appraoch
%  Input Parameters:-------------------------------------------------------
% E:	Imposed macroscale strain rate tensor, a 3-by-3 matrix.
% n:	The number of rheologically distinct elements (RDEs).
% a:	The initial three semi-axes of clasts (a1>a2>a3; a3=1), a 3-by-n matrix.
% q:	The initial orientations of clasts as transformation matrix, 
%       a 3-by-3-by-n matrix. 
% eta:	The effective viscosities of clasts defined at the initial state 
%       (the imposed macroscopic strain-rate state), a column vector.
% eta0: The effective viscosities of matrix defined at the initial state
%       (the imposed macroscopic strain-rate state), a scalar.
% Nk:	The stress exponents of clasts, a column vector.
% Nm:	The stress exponents of the matrix, a scalar.
% Output Parameters:-------------------------------------------------------
% S_bar: Bulk stress
% C_bar: Homogenized stiffness tenor
% s:     The partitioned stresses of clasts
% e:     The Partitioned strain rates of clasts
% vis:   The effective viscosities of clasts defined at the partitioned 
%        local states

%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~] = FourIdentity();

%  strain rate invariants at which the viscosities of RDEs are defiend  
   REF  = Inva(E)*ones(1,n);   
%  matrix stiffness
   C0= 2*eta0*Jd;
   M0 = FourTensorInv(C0);

%  clasts' stiffnesses
   Ck = zeros(3,3,3,3,n);
   for j=1:n
       Ck(:,:,:,:,j) = 2*eta(j)*Jd;
   end
%%  preallocate variables:-------------------------------------------------
  Ak = Ck;
  cA = Ck;
  e = zeros(3,3,n); 
  s = e;
  
  parfor r=1:n   

%           %  describe Cm in clasts' coordinates 
%           c_m = Transform(Cm,q(:,:,r));  
%           %  compute the 4th-order Green tensor(T) in RDEs' coordinates            
%           %T = Tfunction_AGLQ(a(:,r), c_bar, p, ww, 2);
%           T = Tfunction(a(:,r), c_m, Alp1, Bet1, ww1, Alp2, Bet2, ww2,...
%                           Alp3, Bet3, ww3, Alp4, Bet4, ww4,...
%                           Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww, 2);
%                           %(2 for incompressible material)
%           %  Eshelby tensor(SS, PI) in RDEs' coordinates
%           SS   = Contract(Contract(Jd,T),c_bar);

          SS=SnPI_vis(a(:,r));
          %  transform the Eshelby tensor(SS, PI) to the global coordinate
          S=Transform(SS,q(:,:,r)');
          
         for jj=1:1000

             % Ak: Global strain rate partitioning tensor
             t1 = Contract(S,Contract(M0,Ck(:,:,:,:,r))); 
             t2 = Jd-S+Nm*t1;
             t3 = FourTensorInv(t2);
             t4 = Jd+(Nm-1)*S;                                                  
             Ak(:,:,:,:,r) = Contract(t3,t4);     
             % partition strain rate from matrix strain rate
             e(:,:,r) = Multiply(Ak(:,:,:,:,r),E); 
             % the second invariant of the strain rate tensor
             epsilonI = (0.5*contract1(e(:,:,r) ,e(:,:,r) ))^0.5;
             % update viscosity 
             zz       = (epsilonI/REF(r))^((1-Nk(r))/Nk(r))*eta(r) ;
             % check strain rate
             alpha    = abs(zz/eta(r) - 1);
             eta(r)   = zz;
             REF(r)   = epsilonI;
             Ck(:,:,:,:,r)  = 2*eta(r)*Jd; 
             if alpha < 0.001 
                 break
             end
             if jj==1000
                    XX=['warning: clast ',num2str(r),' not converged!!'];
                    disp(XX);
             end
         end
         s(:,:,r) = Multiply(Ck(:,:,:,:,r),e(:,:,r));
         cA(:,:,:,:,r)=Contract(Ck(:,:,:,:,r)-C0,Ak(:,:,:,:,r));
  end

  C_bar = C0 + phi*mean(cA,5);
  S_bar = Multiply(C_bar,E);
  vis = eta;
end