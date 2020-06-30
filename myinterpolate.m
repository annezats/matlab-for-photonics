function[ydata]=myinterpolate(Dataset,ycolnum)%assuming xcolnum =1
%Dataset is your struct of imported data from a csv file
%interpolation to do list:
%0. make loop
%---loop thru # of x values
%1. identify empty data points
%---isnan(mynumber) ((y==NaN doesnt work)
%2. identify the two nearest data points on either side
%---keep variables for previous, current, next
%---for the next loop until find y thats not NaN
%3. do the calculation 
%---y = y1 + ((x - x1) / (x2 - x1)) * (y2 - y1) %linear interpolation
%4. stick the number in
%5. make this a function
xdata= Dataset.data(:,1) % := all rows, 1= which column
ydata=Dataset.data(:,ycolnum)
[numrows,numcols]=size(xdata) %need number of rows for indexing
for i=1:1:(numrows-1) 
    %disp(i)
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
disp(ydata)
if isnan(ydata(numrows,1))
    disp("no value for last value btw")
end
%figure
%plot(xdata,ydata)


