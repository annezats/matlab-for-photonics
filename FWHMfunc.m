           function [FWHM,x,in]= FWHMfunc(xdata,ydata,start_in,end_in)
    %if you are calculating for a dip in reflectance, 
    %you have to flip it to min and flip the < > 
    [y,in]=max(ydata(start_in:end_in));
    in=in+start_in;
    x=xdata(in)%resonance freq
    HM=  y/2;  %0.5;
    for in0= flip(1:1:in-1)
        x1=xdata(in0);
        x0=xdata(in0-1);
        y1=ydata(in0);
        y0=ydata(in0-1);
        if y0 <HM & y1>HM 
            slope=(y1-y0)/(x1-x0);
            f0=((HM-y0)/slope)+x0;
            break
        end
        in0=in0-1;
    end
    for in1= in:1:size(ydata)-1
        x2=xdata(in1);
        x3=xdata(in1+1);
        y2=ydata(in1);
        y3=ydata(in1+1);
        if y2 >HM & y3<HM 
            slope=(y3-y2)/(x3-x2);
            f1=((y2-HM)/slope)+x2;
            break
        end
        in1=in1+1;
    end

    
    FWHM=f1-f0;

end