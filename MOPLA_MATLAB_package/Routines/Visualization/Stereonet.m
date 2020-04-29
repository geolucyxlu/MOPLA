function Stereonet(q)

    %  Equal-area projection

    % compute two spherical angles for three axes
      [a1_ang, a2_ang, a3_ang] = ConvertQ2Angs(q);

    % compute r for equal-area projection, lower hemispheres will be plotted
    % a1
      [~,a1]  = find(a1_ang(2,:));  
      r1(a1)  = sqrt(2) * sin(a1_ang(2,a1)./2);
     
    % a2
      [~,a2]  = find(a2_ang(2,:));
      r2(a2)  = sqrt(2) * sin(a2_ang(2,a2)./2);

    % a3       
      [~,a3]  = find(a3_ang(2,:));
      r3(a3)  = sqrt(2) * sin(a3_ang(2,a3)./2);

    % equal-area projections of a1, a2, a3   

    % a1      
      subplot(1,3,1);
      t = 0 : .01 : 2 * pi;
      P = polarplot(t, ones(size(t)));
      set(P, 'Visible', 'off')
      ax = gca;
      ax.ThetaDir = 'clockwise';
      ax.ThetaZeroLocation = 'right';
      ax.ThetaTick = [0 90 180 270];
      hold on
    % plot red dots 
      polarplot(a1_ang(1,a1),r1(a1),'.r')
      hold off
      title('a1')

    % a2   
      subplot(1,3,2);
      P = polarplot(t, ones(size(t)));
      set(P, 'Visible', 'off')
      ax = gca;
      ax.ThetaDir = 'clockwise';
      ax.ThetaZeroLocation = 'right';
      ax.ThetaTick = [0 90 180 270];
      hold on
    % plot red dots 
      polarplot(a2_ang(1,a2),r2(a2),'.r')
      hold off
      title('a2')

    % a3   
      subplot(1,3,3);
      P = polarplot(t, ones(size(t)));
      set(P, 'Visible', 'off')
      ax = gca;
      ax.ThetaDir = 'clockwise';
      ax.ThetaZeroLocation = 'right';
      ax.ThetaTick = [0 90 180 270];
      hold on
    % plot red dots
      polarplot(a3_ang(1,a3),r3(a3),'.r')
      hold off
      title('a3')
  
end
