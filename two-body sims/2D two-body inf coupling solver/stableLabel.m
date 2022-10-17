function appendage= stableLabel(CSVData)       
        if CSVData(end) == 0
            appendage= '_STABLE';
        else
            appendage= '_UNSTABLE';
        end
    %end
end