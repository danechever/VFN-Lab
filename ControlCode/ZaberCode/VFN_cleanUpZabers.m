%% Close the connections and populate the result accordingly

if port.isvalid
    fclose(port);

    if strcmp(port.Status, 'closed')
        delete(port)
        clear port
        fprintf('Port closed\n')
    else
        fprintf('Port did not close correctly\n')    
    end
else 
    fprintf('port was invalid\n')
end
