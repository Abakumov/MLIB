function [ header ] = get_segy_headers( data, info )
%GET_SEGY_HEADERS get headers from segy structure
%   This script search for available headers and save them in compact
%   structure. Also provides information about available headers. Optional
%   parameter _info_ provides information about data headers. 
%   Abakumov Ivan
%   31st March, 2015
%   Hamburg University
%   e-mail: abakumov_ivan@mail.ru
    if nargin < 1
        error('get_segy_headers: seismic data is required')
    end
    if nargin < 2
        info = 0;     
    end

    [nheaders, ~] = size(data.headers); 
    if info~=0
        disp(data.header_info);
    end
    
    for i=1:nheaders
        % cdp
        if strcmp( char(data.header_info(i,1)), 'cdp')
            header.cdp = data.headers(i, :); 
        end
        % cdpx
        if strcmp( char(data.header_info(i,1)), 'cdp_x')
            header.cdpx = data.headers(i, :); 
        end
        % cdpy
        if strcmp( char(data.header_info(i,1)), 'cdp_y')
            header.cdpy = data.headers(i, :); 
        end
        % offset
        if strcmp( char(data.header_info(i,1)), 'offset')
            header.offset = data.headers(i, :); 
        end
        % sx
        if strcmp( char(data.header_info(i,1)), 'sou_x')
            header.sx = data.headers(i, :); 
        end
        % sy
        if strcmp( char(data.header_info(i,1)), 'sou_y')
            header.sy = data.headers(i, :); 
        end
        % sz
        if strcmp( char(data.header_info(i,1)), 'sou_elev')
            header.sz = data.headers(i, :); 
        end
        % rx
        if strcmp( char(data.header_info(i,1)), 'rec_x')
            header.rx = data.headers(i, :); 
        end
        % ry
        if strcmp( char(data.header_info(i,1)), 'rec_y')
            header.ry = data.headers(i, :); 
        end
        % rz
        if strcmp( char(data.header_info(i,1)), 'rec_elev')
            header.rz = data.headers(i, :); 
        end
         % inline
        if strcmp( char(data.header_info(i,1)), 'iline_no')
            header.inline = data.headers(i, :); 
        end
         % xline
        if strcmp( char(data.header_info(i,1)), 'xline_no')
            header.xline = data.headers(i, :);
        end
        % ffid
        if strcmp( char(data.header_info(i,1)), 'ffid')
            header.ffid = data.headers(i, :); 
        end
        % tracf
        if strcmp( char(data.header_info(i,1)), 'o_trace_no')
            header.tracf = data.headers(i, :); 
        end 
    end
end

