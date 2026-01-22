function navigateNucleiPanels(fig, direction)
% Simulates scrolling by changing the RowHeight of the nuclei main content grid.

handles = guidata(fig);

% Total number of panels that can be scrolled
numPanels = 6; % Spacing, Deconv, Preproc, BG, Segmentation, Filter

% Current state is the index of the top-most visible panel
if ~isfield(handles.lastUsed, 'nucNavPanelIndex')
    handles.lastUsed.nucNavPanelIndex = 0; % Start at Stage 0
end
currentState = handles.lastUsed.nucNavPanelIndex;
newState = currentState;

if strcmp(direction, 'down')
    % Allow 4 states for nuclei (0, 1, 2, 3)
    if newState < 3
        newState = newState + 1;
    end
elseif strcmp(direction, 'up')
    % Stop scrolling up if at the top
    if newState > 0
        newState = newState - 1;
    end
end

if newState == currentState
    return; % No change needed
end

% Get the handle to the main grid for nuclei
mainContentGrid = handles.nucMainContentGrid;

% Define the new row heights based on the state
newRowHeights = repmat({'fit'}, 1, 7); % 6 panels + 1 spacer
newRowHeights{end} = '1x'; % Spacer is always 1x

% Set panel visibility based on navigation state
% Grid rows: 1=Spacing, 2=Deconv, 3=Preproc, 4=BG, 5=Segmentation, 6=Filter, 7=Spacer
if newState == 0
    % Stage 0: Show panels 1-5 (Spacing, Deconv, Preproc, BG, Segmentation)
    visiblePanels = 1:5;
elseif newState == 1
    % Stage 1: Show panels 2-6 (Deconv, Preproc, BG, Segmentation, Filter)
    visiblePanels = 2:6;
elseif newState == 2
    % Stage 2: Show panels 3-6 (Preproc, BG, Segmentation, Filter)
    visiblePanels = 3:6;
elseif newState == 3
    % Stage 3: Show panels 4-6 (BG, Segmentation, Filter)
    visiblePanels = 4:6;
else
    % Fallback: Stage 0
    visiblePanels = 1:5;
end

% Apply visibility to all panels
for i = 1:numPanels
    if ismember(i, visiblePanels)
        newRowHeights{i} = 'fit';
    else
        newRowHeights{i} = 0;
    end
end


mainContentGrid.RowHeight = newRowHeights;

% Update the stored state
handles.lastUsed.nucNavPanelIndex = newState;

% Save the updated handles and parameters (same as channels)
guidata(fig, handles);
% Navigation state is not auto-saved - user must manually export parameters
end
