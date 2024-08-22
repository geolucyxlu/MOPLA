 function [S_bar,C_bar,e,s,vis]= SC_homogenization(E, n, a, q, eta, Ne)
%  SC_homogenization.m       
%  description: To get bulk rehology for polyphase aggragates using
%  Self-consistent homogenization appraoch
%  Input Parameters:-------------------------------------------------------
% E:	Imposed macroscale strain rate tensor, a 3-by-3 matrix.
% n:	The number of rheologically distinct elements (RDEs).
% a:	The initial three semi-axes of RDEs (a1>a2>a3; a3=1), a 3-by-3 matrix.
% q:	The initial orientations of RDEs as transformation matrix, 
%       a 3-by-3-by-n matrix. 
% eta:	The effective viscosities of RDEs defined at the initial state, such
%       as the imposed macroscopic strain-rate state, a column vector.
% Ne:	The stress exponents of RDEs, a column vector.
% Output Parameters:-------------------------------------------------------
% S_bar: Bulk stress
% C_bar: Homogenized stiffness tenor
% e:     Partitioned strain rates of RDEs
% s:     partitioned stresses of RDEs
% vis:   The effective viscosities of RDEs defined at the partitioned local
%        states
%% ------------------------------------------------------------------------ 
% tic
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~] = FourIdentity();
   
%  obtain weights and nodes before the loop 
   gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 

%  Generate the compliance tensors(m), the tangent compliances(mtan)
%  and the pre-strain-rate term(e0)of RDEs
   m      = zeros(3,3,3,3,n);
   mtan   = zeros(3,3,3,3,n);
   e0     = zeros(3,3,n);
   for i= 1:n 
       m(:,:,:,:,i)   = 1/(2*eta(i))*Jd;
       mtan(:,:,:,:,i)= Ne(i)* m(:,:,:,:,i);
       e0(:,:,i)      = (1-Ne(i))*E;
   end
%  strain rate invariants at which the viscosities of RDEs are defiend  
   REF  = Inva(E)*ones(n,1);                           
   
%%  initial guess of the properties of the HEM:----------------------------

%  initial homogenized compliance(M_bar) of HEM
   M_bar = 1/n*sum(mtan,5);                 
%  initial homogenized stiffness(C_bar) of HEM
   C_bar = FourTensorInv(M_bar);
%  initial pre-strain rate term(E0) of HEM   
   mm = FourTensorInv(1/n*sum(m,5));     
   t1 = Contract(M_bar,mm);              
   E0 = Multiply((Jd-t1),E); 
%  macroscale deviatoric stress due to E
   S  = Multiply(C_bar,E-E0);

%  remote deviatoric stress (average deviatoric stress in matrix) due to Em
   Sm = S;        
%  invariant of macroscale deviatoric stress            
   SI = Inva(S);   
%  partitioned deviatoric stress fields inside RDEs
   s  = repmat(S,1,1,n);       
   
%%  preallocate variables:--------------------------------------------------
    e  = zeros(3,3,n);     
    MB   = zeros(3,3,3,3,n);   % for <m:B>
    B    = zeros(3,3,3,3,n);   % for <B>
    Mb   = zeros(3,3,n);       % for <m:beta+e0>
    b    = zeros(3,3,n);       % for <beta>
           
%%  A Self-consistent approach---------------------------------------------

     for i=1:100
            %   Outer loop:When all elements are calculated, using the new set
            %   of M(k),B(k),e0(k),beta(k) update new M_bar and E0 until it
            %   approaches balance.   

    %% parallel computing: Partitioning
       parfor r=1:n
            %   S PI and H_arc for each RDP are calculated starting with
            %   the M_bar and E0 defined at the homogeneous macroscopic
            %   strain-rate state.

            %  describe C in RDEs' coordinates 
            c_bar = Transform(C_bar,q(:,:,r));  
            %  compute the 4th-order Green tensor(T) in RDEs' coordinates            
            %T = Tfunction_AGLQ(a(:,r), c_bar, p, ww, 2);
            T = Tfunction(a(:,r), c_bar, Alp1, Bet1, ww1, Alp2, Bet2, ww2,...
                          Alp3, Bet3, ww3, Alp4, Bet4, ww4,...
                          Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww, 2);
                          %(2 for incompressible material)

            %  Eshelby tensor(SS, PI) in RDEs' coordinates
            SS   = Contract(Contract(Jd,T),c_bar);
            
            %  transform the Eshelby tensor(SS, PI) to the global coordinate
            %  for calculate H_arc and <PI:S^-1>
            ss      = Transform(SS,q(:,:,r)');
            invs    = FourTensorInv(ss);
            %  the 4th-order Interaction tensor(H_arc) 
            H_arc  = Contract(FourTensorInv(invs-Jd),M_bar);   

            for k=1:1000
                %   Inner loop: B and beta(b) for each RDP are calculated 
                %   starting with the M_bar and E0 defined at the macrosclae
                %   strain-rate state uitil is appraches the balance.

                %  the 4th-order stress-partitioning tensor(B)                               
                t1             = FourTensorInv((mtan(:,:,:,:,r) + H_arc));                   
                B(:,:,:,:,r)   = Contract(t1,(M_bar + H_arc));           
                t2             = E0-e0(:,:,r);
                %  the second order stress-partitioning tensor(beta)
                b(:,:,r)       = Multiply(t1,t2);   
                %  calculate new partitioned stress field inside RDE
                sNew           = Multiply(B(:,:,:,:,r),Sm)+ b(:,:,r);                          
                %  compare the current partitioned stress and the last one                 
                alpha          = abs(Inva(sNew-s(:,:,r))/Inva(s(:,:,r)));
                %  replace the last partitioned stress with the current one
                s(:,:,r)       = sNew;    
                %  calculate new partitioned strain rate of RDE               
                e(:,:,r)       = Multiply(m(:,:,:,:,r),sNew);              
                %  calculate new pre-strain rate term of RDE                
                e0(:,:,r)      = (1-Ne(r))*e(:,:,r);
                %  calculate new strain rate invariant
                eI             = Inva(e(:,:,r));
                %  calculate and update new viscosity of RDE at the new strain rate state               
                t3             = 1/Ne(r)-1;
                eta(r)       =(eI/REF(r))^t3 *eta(r); 
                %  replace the last strain rate invariant with the current one                                      
                REF(r)       = eI;                
                %  calculate new compliance tensor of RDE
                m(:,:,:,:,r)   = 1/(2*eta(r))*Jd; 
                mtan(:,:,:,:,r)= Ne(r) * m(:,:,:,:,r);

                %  The inner loop for an RDE terminates when the current stress 
                %  coincides with the previous one within a specific tolerance
                if alpha<0.001 
                   break
                end   
                if k==1000
                    XX=['warning: RDE ',num2str(r),' not converged!!'];
                    disp(XX);
                end

            end

            MB(:,:,:,:,r) = Contract(mtan(:,:,:,:,r), B(:,:,:,:,r));
            Mb(:,:,r)     = Multiply(mtan(:,:,:,:,r), b(:,:,r)) + e0(:,:,r);

        end
    %% Homogenization

        MB1          = sum(MB,5)/n;
        Mb1          = sum(Mb,3)/n;
        BB           = sum(B,5)/n;
        bb           = sum(b,3)/n;
        %  calculate new homogenized compliance for HEM         
        invB      = FourTensorInv(BB);
        New_M_bar = Contract(MB1,invB);  
        %  compare the current homogenized compliance and the previous one        
        delta = Norm(New_M_bar-M_bar)/Norm(M_bar);
        %  replace the previous homogenized compliance with the current one        
        M_bar = New_M_bar;  
        %  calculate new homogenized stiffness for HEM
        C_bar = FourTensorInv(M_bar);  
        %  calculate new pre-strain rate term for HEM
        E0    = Mb1 - Multiply(M_bar,bb);     
        %  calculate new macroscale deviatoric stress for HEM
        S     = Multiply(C_bar, E - E0);              
        %  calculate new secound invariant of the macroscopic deviatoric stress  
        NewSI = Inva(S);  
        %  compare the current macroscale deviatoric stress and the previous one
        gamma = abs(NewSI/SI-1);
        %  replace the last macroscale deviatoric stress invariant with the current one         
        SI    = NewSI; 
        %  calculate new remote stress
        Sm    = Multiply(invB,(S - bb)); 

        %  The outer loop continues until the current homogenized compliance
        %  and macroscopic deviatoric stress coincide with the previous ones 
        %  respectively within specific tolerances        
        if delta<0.02 && gamma<0.02
%             YY1=['Homogenization Converges at ',num2str(i),' step'];
%             disp(YY1); 
%             YY2=['delta=',num2str(delta),' gamma=',num2str(gamma)];
%             disp(YY2) 
            break    
        end
        if i==100
            YY='warning: homogenization is not converged!!';
            disp(YY);
        end
     end
     vis=eta;
     S_bar= S;
%     toc
end
  
   
 