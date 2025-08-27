function VTKData = VTK_read(filename,ext)
% VTK_read - read from VTK file that contains multiple field data.
%
%   VTKData = VTK_read(filename,ext)
%
%   VTKData is a structure currently containing:
%       1. The vtk file headers
%       2. The grid type
%       3. Point data
%       4. Cell structure
%       5. Cell type
%       6. Size of each cell/field data
%       7. Number of field data
%       8. Each field data name, size, and number type
%   
%   The script is very straight forward and does not have many error
%   catching mechanisms. It relies on the user to identify the structure of
%   the VTK file and adjust accordingly if necessary. However, the bulk of
%   the reading script is written here. The user may need to edit the parts
%   of the script that identifies the position of a data header line where
%   the value is written so that the script properly identifies the amount
%   of data stored.
%
%   Copyright 2022 (c) Zhangxi Feng, University of New Hampshire
%% Read multi-data VTK file
% initialize all variables
VTKData.header = cell(4,1);
VTKData.gridType = [];
VTKData.points = [];
VTKData.cells = [];
VTKData.cellTypes = [];
VTKData.cellDataSize = 0;
VTKData.numFieldData = 0;
VTKData.cellData = {};
VTKData.cellDataName = {};
VTKData.cellDataType = {};
%fid = fopen([filename,ext],'r');
fid = fopen(filename,'r');

% 4 lines of header
VTKData.header{1} = fgets(fid);
if (strcmp(VTKData.header{1}(3:5),'vtk') ~= 1)
    disp('Verify vtk file type and header structure');
    return;
end
VTKData.header{2} = fgets(fid);
VTKData.header{3} = fgets(fid);
VTKData.header{4} = fgets(fid);
VTKData.gridType = VTKData.header{4}(9:end);
% flags for searching a title line to find where the numbers are
flag1 = 0;
flag2 = 0;
flag3 = 0;
flag4 = 0;
while ~feof(fid)
    % read header
    str = fgets(fid);
    if (strcmp(str(1:6),'POINTS') == 1)
        % find the end of the number
        for i = 7:length(str)
            if (str(i) ~= ' ' && flag1 == 0)
                flag1 = 1;
            end
            if (str(i) == ' ' && flag1 == 1)
                flag1 = 0;
                numPoints = str2double(str(7:i-1));
                break;
            end
        end
        % initialize the data structure
        VTKData.points = zeros(numPoints,3);
        for i = 1:numPoints
            str = fgets(fid)
            VTKData.points(i,:) = sscanf(str,'%f',[1 3]);
        end
    elseif (strcmp(str(1:5),'CELLS') == 1)
        % find the end of the number
        for i = 6:length(str)
            if (str(i) ~= ' ' && flag2 == 0)
                flag2 = 1;
            end
            if (str(i) == ' ' && flag2 == 1)
                flag2 = 0;
                numCells = str2double(str(6:i-1));
                totVals = str2double(str(i+3:end));
                break;
            end
        end
        numCols = totVals/numCells;
        VTKData.cells = zeros(numCells,numCols);
        for i = 1:numCells
            str = fgets(fid);
            VTKData.cells(i,:) = sscanf(str,'%f',[1 numCols]);
        end
    elseif (strcmp(str(1:10),'CELL_TYPES') == 1)
        % find the end of the number
        for i = 11:length(str)
            if (str(i) ~= ' ' && flag3 == 0)
                flag3 = 1;
            end
            if (str(i) == ' ' && flag3 == 1)
                flag3 = 0;
                numCellType = str2double(str(11:i-1));
                break;
            elseif (i == length(str))
                numCellType = str2double(str(11:end));
            end
        end
        VTKData.cellTypes = zeros(numCellType,1);
        for i = 1:numCellType
            str = fgets(fid);
            VTKData.cellTypes(i) = sscanf(str,'%f');
        end
    else % CELL_DATA section
        if (strcmp(str(1:9),'CELL_DATA') == 1)
            % find the end of the number
            VTKData.cellDataSize = str2double(str(10:end));
            
            % get number of field data and initialize them
            str = fgets(fid);
            VTKData.numFieldData = str2double(str(16:end));
            VTKData.cellDataName = cell(VTKData.numFieldData,1);
            VTKData.cellData = cell(VTKData.numFieldData,1);
            VTKData.cellDataType = cell(VTKData.numFieldData,1);
            cellDataCount = 0;
        else
            % begin reading field data
            % find out how many columns this data has
            cellDataCount = cellDataCount + 1;
            for i = 1:length(str)
                if (str(i) == ' ')
                    VTKData.cellDataName{cellDataCount} = str(1:i-1);
                    for j = i:length(str)
                        if (str(j) ~= ' ' && flag4 == 0)
                            flag4 = 1;
                        end
                        if (str(j) == ' ' && flag4 == 1)
                            flag4 = 0;
                            curDataCol = str2double(str(i:j));
                            VTKData.cellDataType{cellDataCount} = str(j:end);
                            break;
                        end
                    end
                    break;
                end
            end
            % begin scanning data for each field data
            for i = 1:VTKData.cellDataSize
                str = fgets(fid);
                VTKData.cellData{cellDataCount}(i,:) = sscanf(str,'%f',[1 curDataCol]);
            end
        end % end of cell data section
    end % end of distinguishing different data
end % end of file scanning
fclose(fid);
end