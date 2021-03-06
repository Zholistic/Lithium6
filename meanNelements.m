function [ outputArray ] = meanNelements(inputArray, n)
%Takes in an array 'inputArray' and an integer 'n'. Sums successive n 
%elements in the array and outputs the new array (which will be 1/n the
%size)

newindex = 2;

arrayFinal(1) = inputArray(1);
arrayFinal(2) = mean(inputArray(1:n));

for i=n+1:length(inputArray)
    
    if(mod(i,n) == 0)
     %do averaging
     %i
     newindex = newindex + 1;   
     arrayFinal(newindex) = mean(inputArray(i-(n-1):i));                                    
    end

end

%final points (where length of array not divisible by n)
if(mod(length(inputArray),n) ~= 0)
    newindex = newindex+1;
    arrayFinal(newindex) = mean(inputArray(end - (mod(length(inputArray),n)):end));
end



outputArray = arrayFinal;


end