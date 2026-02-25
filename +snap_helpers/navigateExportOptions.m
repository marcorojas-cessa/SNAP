function navigateExportOptions(fig, direction)
% Navigate staged export options using up/down controls.

handles = guidata(fig);

if ~isfield(handles, 'exportItemChecks') || isempty(handles.exportItemChecks)
    return;
end

if ~isfield(handles, 'numExportItems') || isempty(handles.numExportItems)
    handles.numExportItems = 0;
end
if ~isfield(handles, 'exportVisibleSlots') || isempty(handles.exportVisibleSlots)
    handles.exportVisibleSlots = 2;
end
if ~isfield(handles, 'exportNavState') || isempty(handles.exportNavState)
    handles.exportNavState = 0;
end

numItems = handles.numExportItems;
visibleSlots = max(1, handles.exportVisibleSlots);
numStages = max(1, ceil(max(numItems, 1) / visibleSlots));
currentState = handles.exportNavState;

newState = currentState;
switch lower(direction)
    case 'up'
        if currentState > 0
            newState = currentState - 1;
        end
    case 'down'
        if currentState < (numStages - 1)
            newState = currentState + 1;
        end
    otherwise
        % 'refresh' and unknown directions clamp to valid range.
end

newState = max(0, min(newState, numStages - 1));
handles.exportNavState = newState;

startIdx = newState * visibleSlots + 1;
endIdx = min(numItems, startIdx + visibleSlots - 1);

for i = 1:length(handles.exportItemChecks)
    if i <= numItems && i >= startIdx && i <= endIdx
        handles.exportItemChecks(i).Visible = 'on';
        handles.exportItemChecks(i).Layout.Row = i - startIdx + 1;
    else
        handles.exportItemChecks(i).Visible = 'off';
    end
end

if isfield(handles, 'exportStageLabel') && isvalid(handles.exportStageLabel)
    if numItems < 1
        handles.exportStageLabel.Text = '0-0 / 0';
    else
        handles.exportStageLabel.Text = sprintf('%d-%d / %d', startIdx, endIdx, numItems);
    end
end

if isfield(handles, 'exportUpButton') && isvalid(handles.exportUpButton)
    handles.exportUpButton.Enable = onOff(newState > 0);
end
if isfield(handles, 'exportDownButton') && isvalid(handles.exportDownButton)
    handles.exportDownButton.Enable = onOff(newState < (numStages - 1));
end

guidata(fig, handles);
end

function out = onOff(tf)
if tf
    out = 'on';
else
    out = 'off';
end
end
