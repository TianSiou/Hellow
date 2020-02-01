function [  ] = plot_ew_stepBystep( jj, ew_new, sigma, NO_ew_T, ew_new_1, shift, ew_skew_Cayley )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ( ~isempty(ew_new_1) )
    %ew_new = [ ew_new; ew_new_1 ];
    figure(2)
    plot(real(ew_new_1), imag(ew_new_1),'k*','LineWidth',2);
end
                
flag = 1;
                
switch mod(jj, 4)+1
    case {1,6,11}
        if ( NO_ew_T(jj,1) > 0 && flag == 0 )
            figure(1)
            plot(ew_skew_Cayley(1:NO_ew_T(jj,1),jj),'rx','LineWidth',2);
            plot(real(shift(jj,1)),imag(shift(jj,1)),'k^','LineWidth',2);
        end

        figure(2)
        plot(real(ew_new), imag(ew_new),'m+','LineWidth',2);
        plot(real(sigma),imag(sigma),'kd','LineWidth',2);
    case {2,7,12}
        if ( NO_ew_T(jj,1) > 0 && flag == 0 )
            figure(1)
            plot(ew_skew_Cayley(1:NO_ew_T(jj,1),jj),'g*','LineWidth',2);
            plot(real(shift(jj,1)),imag(shift(jj,1)),'k^','LineWidth',2);
        end 
                    
        figure(2)
        plot(real(ew_new), imag(ew_new),'c*','LineWidth',2);
        plot(real(sigma),imag(sigma),'kd','LineWidth',2);
                        
    case {3,8,13}
        if ( NO_ew_T(jj,1) > 0 && flag == 0 )
            figure(1)
            plot(ew_skew_Cayley(1:NO_ew_T(jj,1),jj),'m+','LineWidth',2);
            plot(real(shift(jj,1)),imag(shift(jj,1)),'k^','LineWidth',2);
        end 
                    
        figure(2)
        plot(real(ew_new), imag(ew_new),'rx','LineWidth',2);
        plot(real(sigma),imag(sigma),'kd','LineWidth',2);
                        
    case {4,9,14}
        if ( NO_ew_T(jj,1) > 0 && flag == 0 )
            figure(1)
            plot(ew_skew_Cayley(1:NO_ew_T(jj,1),jj),'c*','LineWidth',2);
            plot(real(shift(jj,1)),imag(shift(jj,1)),'k^','LineWidth',2);
        end 
                    
        figure(2)
        plot(real(ew_new), imag(ew_new),'g+','LineWidth',2);
        plot(real(sigma),imag(sigma),'kd','LineWidth',2);
                        
    otherwise
        if ( NO_ew_T(jj,1) > 0 && flag == 0 )
            figure(1)
            plot(ew_skew_Cayley(1:NO_ew_T(jj,1),jj),'r+','LineWidth',2);
            plot(real(shift(jj,1)),imag(shift(jj,1)),'k^','LineWidth',2);
        end 
                    
        figure(2)
        plot(real(ew_new), imag(ew_new),'mx','LineWidth',2);
        plot(real(sigma),imag(sigma),'kd','LineWidth',2);
                        
end
                   
if ( NO_ew_T(jj,1) > 0 && flag == 0 )
    figure(1)
    drawnow
    hold on
end
                
figure(2)
drawnow
hold on
                
end

