 function [Sm_evl,S_evl,C_bar_evl,Q_evl,A_evl,Qk_evl,Ak_evl]= MOPLA_secondary(L, n, a, q, eta, Ne, steps, tincr, ns, ak, qk, etak, Nk)
%  MOPLA_secondary.m        
%
% This function allows us to simulate the overall mechanical behavior of 
% the heterogeneous rock mass composed of primary RDEs and to track the 
% shape and orientation evolution of all primary and secondary inclusions. 
% Therefore, we can investigate the fabric development at the observed 
% structural and fabric-element level.
%
% Scripts are based on the self-consistent solution of the partitioning and
% homogenization equations.(Jiang, 2013, 2014, 2016; Qu et al., 2016)
% 
% update: remote stress field  and macroscale stress field  are distinct
%         (Oct,22,2015) by Lucy
% update: add parallel computing (parfor) and reduced some unnecessary
%         variables (May,17,2018)by Lucy
% update: add secondary elements (May,29,2018) by Lucy
% upadte: update the Tfunction 30th July, 2019 by Lucy
% update: remote strain rate and vorticity and macroscopic strain rate and
%         vorticity are distinct (Aug,27,2019) Lucy

% Input parameters---------------------------------------------------------
% L:	Imposed macroscale velocity gradient tensor, a 3-by-3 matrix.
% n: 	The number of primary RDEs. 
% a:	The initial three semi-axes of primary RDEs (a1>a2>a3; a3=1), 
%       a 3-by-3 matrix.
% q:	The initial orientations of primary RDEs, a 3-by-3-by-n matrix. 
% eta:	The effective viscosities of primary RDEs defined at the initial state,
%       such as the imposed macroscopic strain-rate state, a 1-by-n vector.
% Ne:	The stress exponents of primary RDEs, a 1-by-n vector.
% ns:	The number of secondary inclusions in each primary RDE.
% ak:	The initial three semi-axes of secondary inclusions,
%       a 3-by-n-by-ns matrix.
% qk:	The initial orientation of secondary inclusions, 
%       a 3-by-3-by-n-by-ns matrix
% etak:	The effective viscosities of secondary inclusions defined at the 
%       initial state, such as the imposed macroscopic strain-rate state,
%       a ns-by-n matrix.
% Nk:	The stress exponents of secondary inclusions, a ns-by-n matrix.
% steps:	The total computational steps.
% tincr:	The time increment of each computational step; the choice of 
%           tincr must ensure that each computational step represents an 
%           infinitesimal deformation (Jiang, 2007).
% Output variables:--------------------------------------------------------
% C_bar_evl:	The evolution with time of the homogeneous macroscopic 
%               stiffness tensor of the heterogeneous rock mass,
%               a 3-by-3-by-3-by-3-by-steps matrix.
% S_bar_evl:	The evolution with time of the macroscopic deviatoric stress
%               tensor due to the imposed macroscopic strain rate,
%               a 3-by-3-by-steps matrix.
% Q _evl:	The orientation evolution of all primary RDEs with time, 
%           a 3-by-3-by-n-by-steps matrix.
% A_evl:	The shape evolution of all primary RDEs with time, 
%           a 3-by-n-by-steps matrix.
% Qk _evl:	The orientation evolution of all secondary inclusions with time,
%           a 3-by-3-by-n-by-steps matrix.
% Ak_evl:	The shape evolution of all secondary inclusions with time, 
%            a 3-by-n-by-steps matrix.
%
%% ------------------------------------------------------------------------
%  Decomposition of the imposed flow L into the strain rate tensor E and the
%  vorticity tensor W        
   E = 0.5 * (L + L');
   W = 0.5 * (L - L');   
%  Generate 4th-order identity tensors   
   [Jd, ~, Ja, ~] = FourIdentity();
   
%  Obtain weights and nodes before the loop 
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
       m(:,:,:,:,i)   = 1/(2*eta(1,i))*Jd;
       mtan(:,:,:,:,i)= Ne(1,i)* m(:,:,:,:,i);
       e0(:,:,i)      = (1-Ne(1,i))*E;
   end
%  strain rate invariants at which the viscosities of primary and secondary
%  inclusions are defiend    
   REF   = Inva(E)*ones(1,n);                           
   REFk  = Inva(E)*ones(ns,n);
   
%%  initial guess of the properties of the HEM:-----------------------------

%  initial homogenized compliance(M_bar) of HEM
   M_bar = 1/n*sum(mtan,5);                 
%  initial homogenized stiffness(C_bar) of HEM
   C_bar = FourTensorInv(M_bar);
%  initial pre-strain rate term(E0) of HEM   
   mm    = FourTensorInv(1/n*sum(m,5));     
   t1    = Contract(M_bar,mm);              
   E0    = Multiply((Jd-t1),E); 
%  the macroscopic deviatoric stress due to E
   S     = Multiply(C_bar,E-E0);     
%  remote deviatoric stress (average strain rate in matrix) 
   Em    = E;
%  remote vorticity (average vorticity in matrix)
   Wm    = W;
%  remote deviatoric stress (average deviatoric stress in matrix) due to Em
   Sm    = S;        
%  invariant of the HEM deviatoric stress(macroscale)             
   SI    = Inva(S);   
%  partitioned stress fields inside RDEs
   s     = repmat(S,1,1,n);       
   
%%  preallocate variables:--------------------------------------------------
        SS = zeros(3,3,3,3,n); 
        PI = zeros(3,3,3,3,n);  
        e  = zeros(3,3,n);     
         
        Q_evl     = zeros(3,3,n,steps);
        A_evl     = zeros(3,n,steps);
        eta_evl   = zeros(1,n,steps);
        C_bar_evl = zeros(3,3,3,3,steps);
        M_bar_evl = zeros(3,3,3,3,steps); 
        S_evl = zeros(3,3,steps);
        Sm_evl = zeros(3,3,steps);
        
        PInS = zeros(3,3,3,3,n);   % for <PI:S^-1>
        MB   = zeros(3,3,3,3,n);   % for <m:B>
        B    = zeros(3,3,3,3,n);   % for <B>
        Mb   = zeros(3,3,n);       % for <m:beta+e0>
        b    = zeros(3,3,n);       % for <beta>
         
        Qk_evl     = zeros(3,3,n,ns,steps);
        Ak_evl     = zeros(3,n,ns,steps);
 %%  A Self-consistent approach---------------------------------------------
 for h=1:steps   
     tic   
     for i=1:100    
            %   Outer loop:When all elements are calculated, using the new set
            %   of M(k),B(k),e0(k),beta(k) update new M_bar and E0 until it
            %   approaches balance.   
            
%% parallel computing: Partitioning
       parfor r=1:n
            %   S PI and H_arc for each RDP are calculated starting with
            %   the M_bar and E0 defined at the homogeneous macrosclae
            %   strain-rate state.

            %  describe C in RDEs' coordinates 
            c_bar = Transform(C_bar,q(:,:,r));            
            %  compute the 4th-order Green tensor(T) in RDEs' coordinates           
            %T = Tfunction_AGLQ(a(:,r), c_bar, p, ww, 2);
            T = Tfunction(a(:,r), c_bar, Alp1, Bet1, ww1,Alp2, Bet2, ww2,...
                             Alp3, Bet3, ww3, Alp4, Bet4, ww4,...
                            Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww, 2);
                            %(2 for incompressible)
            %  Eshelby tensor(SS, PI) in RDEs' coordinates
            SS(:,:,:,:,r)    = Contract(Contract(Jd,T),c_bar);
            PI(:,:,:,:,r)    = Contract(Contract(Ja,T),c_bar);
            %  transform the Eshelby tensor(SS, PI) to the global coordinate
            %  for calculate H_arc and <PI:S^-1>  
            ss      = Transform(SS(:,:,:,:,r),q(:,:,r)');
            pi      = Transform(PI(:,:,:,:,r),q(:,:,r)');
            invs    = FourTensorInv(ss);
            PInS(:,:,:,:,r)  = Contract(pi,invs); 
            %  the 4th-order Interaction tensor(H_arc) 
            H_arc  = Contract(FourTensorInv(invs-Jd),M_bar);    

            for k=1:100
                %   Inner loop: B and beta(b) for each RDP are calculated starting with
                %   the M_bar and E0 defined at the homogeneous macrosclae
                %   strain-rate state uitil is approaches the balance.

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
                e0(:,:,r)      = (1-Ne(1,r))*e(:,:,r);
                %  calculate new strain rate invariant
                eI             = Inva(e(:,:,r));
                %  calculate and update new viscosity of RDE at the new strain rate state               
                t3             = 1/Ne(1,r)-1;
                eta(1,r)       =(eI/REF(:,r))^t3 *eta(1,r); 
                %  replace the last strain rate invariant with the current one                                      
                REF(:,r)       = eI;                
                %  calculate new compliance tensor of RDE
                m(:,:,:,:,r)   = 1/(2*eta(1,r))*Jd; 
                mtan(:,:,:,:,r)= Ne(1,r) * m(:,:,:,:,r);
                
                %  The inner loop for an RDE terminates when the current stress 
                %  coincides with the previous one within a specific tolerance
                if alpha<0.01 
                   break
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
        invB         = FourTensorInv(BB);
        New_M_bar    = Contract(MB1,invB);  
        %  compare the current homogenized compliance and the previous one        
        delta        = Norm(New_M_bar-M_bar)/Norm(M_bar);
        %  replace the previous homogenized compliance with the current one        
        M_bar        = New_M_bar;  
        %  calculate new homogenized stiffness for HEM
        C_bar        = FourTensorInv(M_bar);  
        %  calculate new pre-strain rate term for HEM
        E0           = Mb1 - Multiply(M_bar,bb);     
        %  calculate new homogeneous macroscopic deviatoric stress for HEM
        S            = Multiply(C_bar, E - E0);              
        %  calculate new secound invariant of the macroscopic deviatoric stress  
        NewSI        = Inva(Multiply(Jd,S));  
        %  compare the current macroscopic deviatoric stress and the previous one
        gamma        = abs(NewSI/SI-1);
        %  replace the last macroscopic deviatoric stress invariant with the current one         
        SI           = NewSI; 
        %  calculate new remote stress
        Sm           = Multiply(invB,(S - bb)); 
        %  calculate new remote strain-rate
        Em           = Multiply(M_bar,Sm)+E0;
        %  calculate new remote vorticity
        PInS1        = sum(PInS,5)/n; 
        Wm           = W - Multiply(PInS1,(E-Em));
        %  The outer loop continues until the current homogenized compliance
        %  and macroscopic deviatoric stress coincide with the previous ones 
        %  respectively within specific tolerances        
        if delta<0.02 && gamma<0.02
            break
       end
       
     end
     %  write updated eta to eta_evl
     eta_evl(:,:,h) = eta;
     %  write updated C_bar to C_bar_evl
     C_bar_evl(:,:,:,:,h) = C_bar;
     %  write updated M_bar to M_bar_evl
     M_bar_evl(:,:,:,:,h) = M_bar;
     %  write updated S_bar to S_bar_evl
     S_evl(:,:,h)= S;
     Sm_evl(:,:,h)=Sm;
%%  Parallel Computing: Evolution of RDEs   
     parfor j=1:n
            %  describe Em and Wm in RDEs' coordinates 
            E1 = q(:,:,j)* Em *q(:,:,j)'; 
            W1 = q(:,:,j)* Wm *q(:,:,j)';        
            %  strain rate of RDE
            e1    = q(:,:,j)*e(:,:,j)*q(:,:,j)'; %in RDEs' coordinates
            de    = e1 - E1;
            %  vorticity of RDE        
            u1    = Contract(PI(:,:,:,:,j),FourTensorInv(SS(:,:,:,:,j)));
            we    = Multiply(u1, de)+ W1; 
            wEp   = Wd(a(:,j), W1, e1);        
            %  angular velocity of RDE
            theta = we - wEp; 
            %  update Q using 
             qNew = (RodrgRot(-theta * tincr)) * q(:,:,j); 
            %  update a
             aNew = a(:,j).* exp(diag(e1) * tincr); 
            %  make sure that Q and a are in the descending oreder of a(a1>=a2>=a3)          
             qa   = sortrows([qNew aNew],-4);          
            %  Boudinage if the RDE is too elongated or flattened
            % (a1:a3>25 or a2:a3>25)
             if qa(1,4)/qa(3,4)>25
                qa(1,4)=0.5*qa(1,4);
             end

             if qa(2,4)/qa(3,4)>25
                qa(2,4)=0.5*qa(2,4);
             end
             qa = sortrows(qa,-4);       
             a(:,j)   = qa(:,4);
             q(:,:,j) = qa(1:3,1:3);
             %  write updated Q to Q_evl        
             Q_evl(:,:,j,h) = q(:,:,j);
             %  write updated a to A_evl
             A_evl(:,j,h) = a(:,j);    
             
       %% update secondary elements
        % the flow inside the RDE
        le = e1 + we;

        for k=1:ns
            % write the orientation of the secondary element in RDE's
            % coordinate
            q1 = qk(:,:,j,k) * q(:,:,j)';
            % transfer RDE flow into secondary inclusion's coordinate 
            l = q1 * le * q1';
            % decompose the RDE flow l into strain rate d and vorticity w 
            d = 0.5 * (l + l');
            w = 0.5 * (l - l');
            % the viscosity ratio of secondary and primary inclusions
            vr = etak(k,j)/eta(j);
            % the stiffness tensors of secondary and primary inclusions
            ck = 2*Jd*etak(k,j);
            cm = 2*Jd*eta(j);
            % the 4th-order Eshelby tensors(Sk and Pk) for isotropic
            % material
            [Sk,Pk]= SnPI_vis(ak(:,j,k)); 
            % strain-rate, stiffness tensor of secondary inclusion
            [dk,~] = Ed(Ne(j),Nk(k,j),cm,ck,vr,Sk,d,REFk(k,j),Jd,q1);
            % vorticity of secondary 
            u1    = Contract(Pk,FourTensorInv(Sk));
            u2    = dk - d;
            wk    = Multiply(u1, u2)+ w;
            wkp   = Wd(ak(:,j,k), w, dk);
            % angular velocity of secondary element
            thetak = wk - wkp; 
            % update qk
            qk(:,:,j,k) = (RodrgRot(-thetak * tincr)) * qk(:,:,j,k); 
            %  update ak
            ak(:,j,k)= ak(:,j,k).* exp(diag(dk) * tincr); 
            %  make sure that Q and a are in the descending oreder of a(a1>=a2>=a3)          
            qak    = sortrows([qk(:,:,j,k) ak(:,j,k)],-4);                
            ak(:,j,k)   = qak(:,4);
            qk(:,:,j,k) = qak(1:3,1:3);
            %  write updated qk to Qk_evl        
            Qk_evl(:,:,j,k,h) = qk(:,:,j,k);
            %  write updated ak to Ak_evl
            Ak_evl(:,j,k,h) = ak(:,j,k);  
        end
        
     end
     
       XX=['Step ',num2str(h),' completed'];
       disp(XX);
       toc
     
 end 

end
  
   
 