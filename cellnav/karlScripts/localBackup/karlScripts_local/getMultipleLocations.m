function clickedLocations=getMultipleLocations()
% get a bunch of locations from a patch when you click a spot and hit
% 'space'
clickedLocations = [];

% Enable data cursor mode

%set(dcm, 'UpdateFcn', @cursorUpdateFcn, 'DisplayStyle', 'datatip');

% Set up a key press event handler
set(gcf, 'KeyPressFcn', @keyPressFcn);

% Data cursor update function
% function output_txt = cursorUpdateFcn(~, event_obj)
%     pos = event_obj.Position;
%     output_txt = {['X: ', num2str(pos(1))], ...
%                   ['Y: ', num2str(pos(2))], ...
%                   ['Z: ', num2str(pos(3))]};    
% end

% Key press event handler
function keyPressFcn(~, event)
    if strcmp(event.Key, 'space')
        % Get the current cursor info
        cursorInfo = getCursorInfo(dcm);
        if ~isempty(cursorInfo)
            % Extract the clicked location and add it to the array
            clickedLocation = cursorInfo(1).Position;
            clickedLocations = [clickedLocations; clickedLocation];
            fprintf('Clicked location added: [X: %.2f, Y: %.2f, Z: %.2f]\n', ...
                clickedLocation(1), clickedLocation(2), clickedLocation(3));
        end
    end
end
end