function Out = aacgm(lat,lon,time,alt,varargin)
% MATLAB wrapper for the AACGM model and MEX function
% To compile MEX file simply call this function
% By default this function converts from geodetic to aacgm coordinates.
% To convert from magnetic to geodetic, use the optional 'mag2geo' keyword.
%
% INPUTS
%   latitude(deg), longitude(deg), time(datenum or datetime), altitude(km)
%
% OUTPUTS
%   structure with fields
%   the "error" field is a 2xN matrix of the return codes from AACGM:
%   1xN: AACGM_v2_SetDateTime
%   2xN: AACGM_v2_Convert
%
% time, lon, lat, alt, mlt(hours), mlat(deg), mlon(deg), malt (km)
% OPTIONAL INPUTS
%
%   'geo2mag' or 'G2A': convert from geodetic to aacgm (default)
%
%   'mag2geo' or 'A2G': convert from aacgm to geodetic
%
%   'sat' mode: in this mode, all inputs must be the same size, including
%   altitude
%
%   'map' mode: in this mode, all unique values of inputs are expanded in a
%   4d grid. The output Ne is of a size [n lat,n lon,n alt,n time]
%
%   'TRACE', 'ALLOWTRACE': by default aacgm is only defined up to 2000km.
%    these options allow field line tracing to extend above this limit, but
%    is very slow. 'TRACE' forces field tracing, whereas ALLOWTRACE only
%    enables tracing if needed.
%
%   Ben Reid 2020, 2023
%   Modifed for E-CHAIM - Anthony McCaffrey 2023

%Directory which contains the mex file,
%the "source" directory, and the "coeffs" directory
RootDir = fileparts(which("aacgm"));

%get the relavent directories
aacgm_dirs = init(fullfile(RootDir,"coeffs"));
ModelDir = fullfile(RootDir,'source');

ErrorStr = 'MATLAB:AACGM';


if ~exist(fullfile(RootDir,'mex_aacgm'),'file')
    %COMPFLAGS='$COMPFLAGS -D__USE_MINGW_ANSI_STDIO'
    mex('-v',...
        ['-I',ModelDir],...
        '-output',fullfile(RootDir,'mex_aacgm'), ...
        fullfile(ModelDir,'aacgmlib_v2.c'), ...
        fullfile(ModelDir,'mlt_v2.c'), ...
        fullfile(ModelDir,'igrflib.c'), ...
        fullfile(ModelDir,'rtime.c'), ...
        fullfile(ModelDir,'astalglib.c'), ...
        fullfile(RootDir,'mex_aacgm.c'));
     
end

code = 0; % default is G2A
mode='default';
i=1;
while i <= numel(varargin)
    if ischar(varargin{i}) || isstring(varargin{i})
        if strcmpi(varargin{i},'default')
            mode = default;
        elseif strcmpi(varargin{i},'sat')
            mode = 'sat';
        elseif strcmpi(varargin{i},'map')
            mode='map';
        elseif strcmpi(varargin{i},'geo2mag') || strcmpi(varargin{i},'G2A')
            code = code + 0;
        elseif strcmpi(varargin{i},'mag2geo') || strcmpi(varargin{i},'A2G')
            code = code + 1;
        elseif strcmpi(varargin{i},'TRACE')
            code = code + 2; % /* use field-line tracing to compute coordinates     */
        elseif strcmpi(varargin{i},'ALLOWTRACE')
            code = code + 4;  % /* if height is >2000 km use tracing, else use coefs */
        elseif strcmpi(varargin{i},'BADIDEA')
            code = code + 8;     % /* use coefficients above 2000 km; Terrible idea!!   */
        else
            error(ErrorStr,[...
                'Invalid option: ',varargin{i}]);
        end
    else
        error(ErrorStr,[...
            'Invalid option.'])
    end
    i=i+1;
end


if any((lat(:)<-90) | (lat(:)>90))
    error(ErrorStr,[...
        'Input latitude out of range [-90,90]'])
end

if any((lon(:)<-180) | (lon(:)>360))
    error(ErrorStr,[...
        'Input longitude out of range [-180,360]'])
end

time = datetime(datenum(time),'ConvertFrom','datenum');

if strcmpi(mode,'default')
    if (numel(time) ~=1) && (numel(time) ~= numel(lat))
        error(ErrorStr,[...
            'In "default" mode, the number of input times must be either',...
            ' equal to one, or equal to the number of lat/lon pairs'])
    end

    tvec=[numel(lat),numel(lon)];
    tvec=tvec==tvec';
    if ~all(tvec(:))
        error(ErrorStr,[...
            'In "default" mode lat and lon must have the same size'])
    end
    if numel(time) == 1
        time = repmat(time, size(lon));
    end
    Out.time = time(:);
    Out.lon = lon(:);
    Out.lat = lat(:);
    Out.alt = alt(:);

    OutFmt = @(X) reshape(X,numel(lon),numel(alt));

    [lon,~] = ndgrid(lon(:),alt(:));
    [lat,~] = ndgrid(lat(:),alt(:));
    [time,alt] = ndgrid(time(:),alt(:));


elseif strcmpi(mode,'sat')
    tvec=[numel(lat),numel(lon),numel(alt),numel(time)];
    tvec=tvec==tvec';
    if ~all(tvec(:))
        error(ErrorStr,[...
            'In "sat" mode lat, lon, alt, time must all have the same size'])
    end

    Out.time = time(:);
    Out.lon = lon(:);
    Out.lat = lat(:);
    Out.alt = alt(:);
    OutFmt = @(X) X(:);

elseif strcmpi(mode,'map')

    lat = unique(lat);
    lon = unique(lon);
    alt = unique(alt);
    time = unique(time);

    Out.lat = lat(:);
    Out.lon = lon(:);
    Out.alt = alt(:);
    Out.dates = time(:);
    OutFmt = @(X) reshape(X,numel(lat),numel(lon),numel(alt),numel(time));
    [lat,lon,alt,time] = ndgrid(lat(:),lon(:), alt(:), time(:));
end

[d, Out.mlt, err] = mex_aacgm(char(aacgm_dirs(1)), char(aacgm_dirs(2)), ...
    double([lat(:),lon(:),alt(:)])', ...
    datevec(time)',double(code));

if bitand(code,1)
    Out.mlat = Out.lat;
    Out.mlon = Out.lon;
    Out.malt = Out.alt;
    Out.lat = OutFmt(d(1,:))';
    Out.lon = OutFmt(d(2,:))';
    Out.alt = OutFmt(d(3,:))';
else
    Out.mlat = OutFmt(d(1,:))';
    Out.mlon = OutFmt(d(2,:))';
    Out.malt = OutFmt(d(3,:))';
end

Out.mlt = OutFmt(Out.mlt)';

Out.errors = err;

end

function vars = init(coefs_dir)
    
    %get the list of files
    f = dir(coefs_dir);

    %get the desired files
    names = arrayfun(@(x) string(f(x).name),1:length(f),'un',1);
    ind(1) = find(contains(names,"aacgm_coeffs-"),1,'first');
    ind(2) = find(contains(names,"igrf"),1,'first');

    %check that the files were found
    if isempty(ind(1)); warning("AACGM Coefficients files not found."); end
    if isempty(ind(2)); warning("IGRF Coefficients files not found."); end

    %split the prefix from the full file name
    t = strsplit(names(ind(1)),'-');
    t(3) = "";
    prefix = strjoin(t,'-');

    %set the env vars
    vars(1) = fullfile(coefs_dir,prefix);
    vars(2) = fullfile(coefs_dir,names(ind(2)));
    setenv("AACGM_v2_DAT_PREFIX",vars(1));
    setenv("IGRF_COEFFS",vars(2));

end
