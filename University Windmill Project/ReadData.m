% determines useable data from files and stores into an array


function [outputArray] = ReadData(fileName, dim)
    data = readmatrix(fileName);
    if dim == 1
        outputArray = [];
        counter = 1;
        for j = 1:length(data(1,:))
            for i = 1:length(data)
                % if we have a useful value, add to the 1D array.
                if ~isnan(data(i,j))
                    outputArray(counter) = data(i,j);
                    counter = counter+1;
                end
            end
        end
    else
        outputArray = readtable(fileName);
    end
end