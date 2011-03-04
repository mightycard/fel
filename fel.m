%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    irpss=load(['data/tau_',int2str(i),'.dat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    x=irpss(:,1);
    z=irpss(:,3);

    ex=irpss(:,4);
    ex1=irpss(:,9);
    ex2=irpss(:,10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%    subplot(2,2,1);
%    ititle=(['E_x ',int2str(i)]);
%    plot3(x/pi,z/(2.0*pi),ex2),title(ititle);
%    xlabel('pi*x/a (mm)'),ylabel('z/k_z '),zlabel('Ex (V/m)');
%    axis([0 1 0 2 -1 1]);
    
%    subplot(2,2,3);
%    ititle=(['E_{0x} ',int2str(i)]);
%    plot3(x/pi,z/(2.0*pi),ex2),title(ititle);
%    xlabel('pi*x/a (mm)'),ylabel('z/k_z '),zlabel('E0x (V/m)');
%    axis([0 1 0 2 -1 1]);
    
%    subplot(2,2,4);
    ititle=(['E_{1x} ',int2str(i)]);
    plot3(z/(2.0*pi),x/pi,ex2),title(ititle);
    xlabel('z/k_z '),ylabel('pi*x/a (mm)'),zlabel('E1x (V/m)');
%    axis([0 1 0 2 -1 1]);
       
    filename=(['img/ex',int2str(i),'.jpg']);
    saveas(gcf,filename,'jpg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%