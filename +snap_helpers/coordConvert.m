classdef coordConvert
% COORDCONVERT - Coordinate conversion utilities for SNAP
%
% SNAP uses ARRAY CONVENTION for storing coordinates internally:
%   [row, col, slice] where:
%   - row: 1st dimension (vertical index, Y-axis in Cartesian)
%   - col: 2nd dimension (horizontal index, X-axis in Cartesian)
%   - slice: 3rd dimension (depth index, Z-axis in Cartesian)
%
% For DISPLAY (plots, graphs), we convert to CARTESIAN CONVENTION:
%   [x, y, z] where:
%   - x: horizontal (col)
%   - y: vertical (row)
%   - z: depth (slice)
%
% Array access: imageData(row, col, slice) - DIRECT, no conversion
% Display plots: use [x, y, z] = [col, row, slice]

    methods (Static)
        function [x, y, z] = array2cartesian(row, col, slice)
            % Convert ARRAY convention to CARTESIAN for display
            % 
            % Inputs:
            %   row   - 1st dimension (vertical index)
            %   col   - 2nd dimension (horizontal index)
            %   slice - 3rd dimension (depth index)
            %
            % Outputs:
            %   x - horizontal position (for plots)
            %   y - vertical position (for plots)
            %   z - depth position (for plots)
            
            x = col;
            y = row;
            z = slice;
        end
        
        function [row, col, slice] = cartesian2array(x, y, z)
            % Convert CARTESIAN (display) to ARRAY convention for storage
            %
            % Inputs:
            %   x - horizontal position
            %   y - vertical position
            %   z - depth position
            %
            % Outputs:
            %   row   - 1st dimension (for array indexing)
            %   col   - 2nd dimension (for array indexing)
            %   slice - 3rd dimension (for array indexing)
            
            row = y;
            col = x;
            slice = z;
        end
        
        function cart_coords = array2cartesian_vec(array_coords)
            % Convert ARRAY convention coordinates to CARTESIAN (vectorized)
            %
            % Input:
            %   array_coords - Nx3 matrix of [row, col, slice]
            %
            % Output:
            %   cart_coords - Nx3 matrix of [x, y, z]
            
            cart_coords = array_coords(:, [2, 1, 3]);  % Swap col and row
        end
        
        function array_coords = cartesian2array_vec(cart_coords)
            % Convert CARTESIAN coordinates to ARRAY convention (vectorized)
            %
            % Input:
            %   cart_coords - Nx3 matrix of [x, y, z]
            %
            % Output:
            %   array_coords - Nx3 matrix of [row, col, slice]
            
            array_coords = cart_coords(:, [2, 1, 3]);  % Swap x and y
        end
    end
end

