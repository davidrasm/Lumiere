function [ output_args ] = batch_treeSim()
% Simulate a batch of trees under same params

start = 1;
sims = 2;

file_prefix = 'fitBD_randParams_2sites'; % to name files
xml_template = 'fitBD_randParams_2sites_template.xml'; % xml template

for s = start:sims
    
    s
    
    file_out = strcat(file_prefix,'_sim',num2str(s)); %file_prefix;
    
    %write_params_file(file_out); % Moved this inside sim while loop so don't get stuck
    
    params_file = strcat(file_out,'_params');
    main_fitBD_sim(params_file, file_out);
    
    convertPhyForBDMM(file_out);
    
    writeBeastXMLRandomParams(file_out, file_out, xml_template, params_file);
    
end


end

