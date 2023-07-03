function gif2mp4(gifFilePath,varargin)


%%% how to use it

% gifFilePath = 'input.gif';
% gif2mp4(gifFilePath);  % Default frame rate (10 fps) and same file name for output
% gif2mp4(gifFilePath, 'FrameRate', 15);  % Custom frame rate (15 fps) and same file name for output
% gif2mp4(gifFilePath, 'output.mp4');  % Default frame rate (10 fps) and custom output file name
% gif2mp4(gifFilePath, 'output.mp4', 'FrameRate', 24);  % Custom frame rate (24 fps) and custom output file name



% Check if the second argument is provided
if nargin < 2 
    % Set mp4FilePath equal to gifFilePath
    mp4FilePath = gifFilePath;
    [~, fileName, ~] = fileparts(mp4FilePath);
    mp4FilePath = [fileName '.mp4'];  % Append '.mp4' extension
end

% Check if the frame rate option is provided
frameRate = 10;  % Default frame rate
if nargin > 2 && rem(numel(varargin), 2) == 0
    mp4FilePath = gifFilePath;
    [~, fileName, ~] = fileparts(mp4FilePath);
    mp4FilePath = [fileName '.mp4']; 
    optionIndex = find(strcmpi(varargin, 'FrameRate'), 1);
    if ~isempty(optionIndex)
        frameRate = varargin{optionIndex + 1};
    end
elseif nargin > 1 && ~rem(numel(varargin), 2) == 0
    mp4FilePath=varargin{1};
    optionIndex = find(strcmpi(varargin, 'FrameRate'), 1);
    if ~isempty(optionIndex)
        frameRate = varargin{optionIndex + 1};
    end
end


% Read the GIF file
[gifFrames, map] = imread(gifFilePath+".gif", 'frames', 'all');
[~, ~, ~, numFrames] = size(gifFrames);
% Create a VideoWriter object with MPEG-4 compression
writerObj = VideoWriter(mp4FilePath, 'MPEG-4');

% Set the frame rate (optional)
writerObj.FrameRate = frameRate;  % Adjust this value as needed

% % Check if frame rate option is provided and update the value if so
% if ~isempty(varargin)
%     numOptions = numel(varargin);
%     if rem(numOptions, 2) ~= 0
%         error('Optional arguments must be provided as name-value pairs.');
%     end
%     for i = 1:2:numOptions
%         optionName = varargin{i};
%         optionValue = varargin{i+1};
%         
%         % Update the value of the frame rate option
%         if strcmpi(optionName, 'FrameRate')
%             writerObj.FrameRate = optionValue;
%         else
%             error('Unknown option: %s', optionName);
%         end
%     end
% end



% Open the VideoWriter object
open(writerObj);

% Write each frame to the MP4 file
for i = 1:1:numFrames
    % Convert indexed image to RGB
    frame = ind2rgb(gifFrames(:, :, i), map);
    
    % Write the frame to the MP4 file
    writeVideo(writerObj, frame);
end

% Close the VideoWriter object
close(writerObj);

% Display a message when the conversion is complete
disp('GIF to MP4 conversion is complete!');


end





