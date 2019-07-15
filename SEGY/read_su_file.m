function data = read_su_file(sufilename)

% Convert SU file to SEGY and read SEGY
% require Seismic Unix package
% Abakumov Ivan
% 25th Feb 2015
% Hamburg University

% if SU file exist, 
if exist(sufilename, 'file') == 2
    % analize SU filename and make proper SEGY filename
    [pathstr,name,ext] = fileparts(sufilename); 
    if strcmp(ext,'.su')
        segyfilename = [pathstr,'/', name, '.sgy']; 
    else
        segyfilename = [pathstr,'/', name, ext, '.sgy']; 
    end

    % if this SEGY already exists:
    if exist(segyfilename, 'file') == 2
        SUFileInfo   = dir(sufilename);
        SEGYFileInfo = dir(segyfilename);
        % if file is new, just read it
        if SEGYFileInfo.datenum > SUFileInfo.datenum
            data = read_segy_file(segyfilename);
        else
            % else delete old SEGY file and create new one
            disp('WARNING: SEGY is old!');
            setenv('SUFNAME', sufilename); 
            setenv('SEGYFNAME', segyfilename); 
            !echo $SUFNAME
            !echo $SEGYFNAME
            ! rm -f $SEGYFNAME
            !segyhdrs < $SUFNAME
            !segywrite < $SUFNAME tape=$SEGYFNAME
            !rm -f header
            !rm -f binary
            % and then read SEGY
            data = read_segy_file(segyfilename);
        end
    else
        % else first convert SU to SEGY
        setenv('SUFNAME', sufilename); 
        setenv('SEGYFNAME', segyfilename); 
        !echo $SUFNAME
        !echo $SEGYFNAME
        !segyhdrs < $SUFNAME
        !segywrite < $SUFNAME tape=$SEGYFNAME
        !rm -f header
        !rm -f binary
        % and then read SEGY
        data = read_segy_file(segyfilename); 
    end
else
    disp('ERROR: SU file does not exist!');
end

