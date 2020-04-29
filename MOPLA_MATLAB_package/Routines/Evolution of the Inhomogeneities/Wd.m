function wp = Wd(a,w,e)
% Wd.m
%
% the Ellipsoid vorticity referred to the frame tracking its semi-axes
%--------------------------------------------------------------------------
 
  wp = zeros(3,3);
  
  for i = 1:3
      for j = 1:3
          if i == j
              wp(i,j) = 0;
          elseif a(i) == a(j)
              wp(i,j) = w(i,j);
          else
              r       = (a(i)^2+a(j)^2)/(a(i)^2-a(j)^2);
              wp(i,j) = r * e(i,j);
          end
      end
  end
  
end
