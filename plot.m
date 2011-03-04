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
    ititle=(['IRPSS (z=0.0',int2str(i),'cm, Q=0.001nC)']);
    subplot(2,2,1);plot(1000*x0,1000*xp0,'r+'),title(ititle);
    axis([-1.5, 1.5, -400, 400]),grid on;
    xlabel('x (mm)'),ylabel('xp (mrad)');
    
    ptitle1=(['PARMELA (z=0.0',int2str(i),'cm, Q=0.001nC)']);
    subplot(2,2,2);plot(10*x1,xp1,'g*'),title(ptitle1);
    axis([-1.5, 1.5, -400, 400]),grid on;
    xlabel('x (mm)'),ylabel('xp (mrad)');
        
    etitle=(['ES (z=0.0',int2str(i),'cm, Q=0.001nC)']);
    subplot(2,2,3);plot(1000*x2,1000*xp2,'k+'),title(etitle);
    axis([-1.5, 1.5, -400, 400]),grid on;
    xlabel('x (mm)'),ylabel('xp (mrad)');
    
    atitle1=(['ASTRA (image=on, bt=on) (z=0.0',int2str(i),'cm, Q=0.001nC)']);
    subplot(2,2,4);plot(1000*x21,1000*px21./pz21,'bx'),title(atitle1);
    axis([-1.5, 1.5, -400, 400]),grid on;
    xlabel('x (mm)'),ylabel('xp (mrad)');
        
    filename=(['img/exvsx',int2str(i),'.jpg']);
    saveas(gcf,filename,'jpg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%