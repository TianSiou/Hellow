
function [] = output(fid, ew_order, EPIteNo, CPUtime, ew, rsdl)

          if ( imag(ew) > 0 )  
             fprintf(fid,'itno(kk,%3.0f) = %4.0f; CPUtime(kk,%3.0f) = %8.1f; ew(kk,%3.0f) = %24.16e+%23.16ei; rsdl(kk,%3.0f) = %13.4e; \n', ...
                     ew_order, EPIteNo, ew_order, CPUtime, ew_order, real(ew), imag(ew), ew_order, rsdl); 
          elseif ( imag(ew) < 0 ) 
             fprintf(fid,'itno(kk,%3.0f) = %4.0f; CPUtime(kk,%3.0f) = %8.1f; ew(kk,%3.0f) = %24.16e-%23.16ei; rsdl(kk,%3.0f) = %13.4e; \n', ...
                     ew_order, EPIteNo, ew_order, CPUtime, ew_order, real(ew), -imag(ew), ew_order, rsdl); 
          else
             fprintf(fid,'itno(kk,%3.0f) = %4.0f; CPUtime(kk,%3.0f) = %8.1f; ew(kk,%3.0f) = %24.16e; rsdl(kk,%3.0f) = %13.4e; \n', ...
                     ew_order, EPIteNo, ew_order, CPUtime, ew_order, ew, ew_order, rsdl); 
          end
 
