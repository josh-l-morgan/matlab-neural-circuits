
swc_root_path = "\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\KxR_P11LGN\CellNav_KxR\Volumes\HighRes2023\Analysis\swc\"
swc_save_path = "D:\swc\"
swc_files = dir(swc_root_path+"*.swc")
for i=1:length(swc_files)
    swc_filename = swc_root_path+swc_files(i).name;
    fid = fopen(swc_filename, 'r');
    neuron = struct('ID', [], 'type', [], 'x', [], 'y', [], 'z', [], 'radius', [], 'parent', []);
    while ~feof(fid)
        % Read a line
        line = fgetl(fid);
        
        % Skip comments and empty lines
        if isempty(line) || line(1) == '#'
            continue;
        end
        
        % Parse the line
        data = sscanf(line, '%d %d %f %f %f %f %d');
        
        % Append the data to the neuron structure
        neuron.ID(end+1) = data(1);
        neuron.type(end+1) = data(2);
        neuron.x(end+1) = data(3);
        neuron.y(end+1) = data(4);
        neuron.z(end+1) = data(5);
        neuron.radius(end+1) = data(6);
        neuron.parent(end+1) = data(7);
    end
    fclose(fid);
    save(swc_save_path+swc_files(i).name(1:end-4)+".n", "neuron", "-v7.3")
%     file_id = fopen(swc_save_path+swc_files(i).name, 'w');
%     fclose(file_id);
end
