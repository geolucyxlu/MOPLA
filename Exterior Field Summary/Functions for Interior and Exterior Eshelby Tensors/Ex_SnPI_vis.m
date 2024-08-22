function [S_Ex,PI_Ex,LAMBDA_Ex,G] = Ex_SnPI_vis(a,ep)
%EX_SNPI_VIS Summary of this function goes here
%   Calculate the exterior-point Eshleby tensors (S, PI, LAMBDA) for
%   isotropic viscous material.
%   The explicit formula for the exterior-point symmetric Eshelby tensor S 
%   (G in Ju and Sun 1999) in isotropic elastic materials has been extended
%   to isotropic viscous materials by Jiang (2016). This MATLAB function
%   uses the equation in Jiang, 2016.
%
%   input: a: semi-axes [a1;a2;a3]
%          ep: exterior points (3*n matrix)


    % number of the exterior points
      n = length(ep);

    % G (Ju and Sun 1999, Eq.4 setting mu=0.5) 
      G = Ex_Gtensor(a,ep);

    % LMABDA_Ex (Jiang 2016 Eq.25 LAMBDA_Ex, S_Ex, and PI_Ex)
      LAMBDA_Ex   = zeros(3,3,n);
      for i=1:3
          for j=1:3
              LAMBDA_Ex(i,j,:) = -1/3.*(G(i,j,1,1,:)+ G(i,j,2,2,:)+ G(i,j,3,3,:));
              LAMBDA_Ex(j,i,:) = LAMBDA_Ex(i,j,:);
          end
      end
      for i=1:3
          t                = 1/3.*squeeze(G(1,1,i,i,:)+ G(2,2,i,i,:)+ G(3,3,i,i,:));
          LAMBDA_Ex(i,i,:) = squeeze(LAMBDA_Ex(i,i,:))+t;
      end

    % S_Ex & PI_Ex
      S_Ex   = zeros(3,3,3,3,n);
      PI_Ex  = zeros(3,3,3,3,n);
      delt   = eye(3);
      for i=1:3
          for j=i:3
              for k=1:3
                  for l=k:3
                      %S_Ex
                      S_Ex(i,j,k,l,:)  = squeeze(G(i,j,k,l,:))+ delt(k,l).*...
                                         squeeze(LAMBDA_Ex(i,j,:));
                      S_Ex(j,i,k,l,:)  = S_Ex(i,j,k,l,:);
                      S_Ex(j,i,l,k,:)  = S_Ex(i,j,k,l,:);
                      S_Ex(i,j,l,k,:)  = S_Ex(i,j,k,l,:);
                      %PI_Ex
                      PI_Ex(i,j,k,l,:) = 1/2.*(delt(j,k).*LAMBDA_Ex(i,l,:) +...
                                         delt(j,l).*LAMBDA_Ex(i,k,:) - delt(i,k)...
                                         .*LAMBDA_Ex(j,l,:) - delt(i,l).*...
                                         LAMBDA_Ex(j,k,:));
                      PI_Ex(i,j,l,k,:) = PI_Ex(i,j,k,l,:);
                      PI_Ex(j,i,k,l,:) = -PI_Ex(i,j,k,l,:);
                      PI_Ex(j,i,l,k,:) = -PI_Ex(i,j,k,l,:);
                  end
              end
          end
      end
      S_Ex(1,1,1,1,:) = -(S_Ex(1,1,2,2,:)+S_Ex(1,1,3,3,:));
      S_Ex(2,2,2,2,:) = -(S_Ex(2,2,1,1,:)+S_Ex(2,2,3,3,:));
      S_Ex(3,3,3,3,:) = -(S_Ex(3,3,1,1,:)+S_Ex(3,3,2,2,:));

end

