function[ActualData]=myinterpolate2(Dataset)
%Dataset is your struct of imported data from a csv file
%^using the importdata('name.csv') function
ActualData=Dataset.data(:,:)
[numrows,numcols]=size(ActualData)
xdata= Dataset.data(:,1) % := all rows, 1= which column
for j=2:1:(numcols)
    ydata=ActualData(:,j)
    for i=1:1:(numrows-1) 
        y=ydata(i,1)
        x=xdata(i,1)
        if isnan(y) %do things if we need to add a value
            if i==1 %edge case for prev values
                x1=0
                y1=0
            else %setting previous values
                x1= xdata((i-1),1)
                y1= ydata((i-1),1)
            end
            n=i+1 %doing things for next values
            y2=ydata(n,1)
            while isnan(y2) %keep going until u get a numerical y
                n=n+1
                y2=ydata(n,1)
            end
            x2=xdata(n,1)
            %time for the formula itself
            ydata(i,1) = y1 + ((x - x1) / (x2 - x1)) * (y2 - y1) %linear interpolation
        end
    end
    ActualData(:,j)=ydata
end
end
