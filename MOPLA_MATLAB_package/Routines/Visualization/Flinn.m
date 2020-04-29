function Flinn(A)

    % Flinn diagram
    
    x = log(A(2,:)./A(3,:));
    y = log(A(1,:)./A(2,:));
    
    
    plot(x,y,'.r',0:8,0:8,'-k')
    xlabel('ln(a2/a3)')
    ylabel('ln(a1/a2)')
    title('Flinn diagram')
    axis square
    
end