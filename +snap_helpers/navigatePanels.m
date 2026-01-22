function navigatePanels(fig, channel_idx, direction)
% Simulates scrolling by changing the RowHeight of the main content grid.

handles = guidata(fig);

% Total number of panels that can be scrolled
numPanels = 8; % Spacing, Deconv, Preproc, BG, Maxima, Gaussian Fit, Fit Filter, Classification

% Current state is the index of the top-most visible panel
currentState = handles.lastUsed.navPanelIndex{channel_idx};
newState = currentState;

if strcmp(direction, 'down')
    % Allow 6 states for channels (0, 1, 2, 3, 4, 5)
    if newState < 5
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

% Get the handle to the main grid for this channel
mainContentGrid = handles.tabMainGrids(channel_idx);

% Define the new row heights based on the state
newRowHeights = repmat({'fit'}, 1, 9); % 8 panels + 1 spacer
newRowHeights{end} = '1x'; % Spacer is always 1x

% Set panel visibility based on navigation state
% Grid rows: 1=Spacing, 2=Deconv, 3=Preproc, 4=BG, 5=Maxima, 6=Gaussian Fit, 7=Fit Filter, 8=Classification, 9=Spacer
if newState == 0
    % Stage 0: Show panels 1-5 (Spacing, Deconv, Preproc, BG, Maxima)
    visiblePanels = 1:5;
elseif newState == 1
    % Stage 1: Show panels 2-6 (Deconv, Preproc, BG, Maxima, Gaussian Fit)
    visiblePanels = 2:6;
elseif newState == 2
    % Stage 2: Show panels 3-7 (Preproc, BG, Maxima, Gaussian Fit, Fit Filter)
    visiblePanels = 3:7;
elseif newState == 3
    % Stage 3: Show panels 4-8 (BG, Maxima, Gaussian Fit, Fit Filter, Classification)
    visiblePanels = 4:8;
elseif newState == 4
    % Stage 4: Show panels 6-8 (Gaussian Fit, Fit Filter, Classification)
    visiblePanels = 6:8;
elseif newState == 5
    % Stage 5: Show panels 7-8 (Fit Filter, Classification) - maximum Classification view
    visiblePanels = 7:8;
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

% Store the new state
handles.lastUsed.navPanelIndex{channel_idx} = newState;

% Save the updated handles
guidata(fig, handles);
% Navigation state is not auto-saved - user must manually export parameters

end
