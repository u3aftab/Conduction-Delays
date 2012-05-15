 function y = in_signal_1d(loc,ampl,sig,nn)
 % !!! REM sig should be smaller than nn/6
 % assigns an input signal with gaussian shape
    if sig>nn/6;
       disp('Warning: The width of input should be less than nn/6 !!!');
    end
    y=zeros(nn,1);
    for i=floor(loc-3*sig):ceil(loc+3*sig);
         imod=mod(i-1,nn)+1;
         dis=(i-loc)^2/(2.*sig^2);
         y(imod)=ampl*exp(-dis);
    end
 return