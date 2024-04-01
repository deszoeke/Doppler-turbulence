function result=nc_varget_lapxm(filename,varname,varargin)
% Read LapXM netcdf files that have missing_value attributes erroneously 
% set to character strings. Convert the character string into the type of
% the variable, and save the string as a new missing_string attribute.
%
% (c) 2009-04-21 Simon de Szoeke

try
    result=nc_varget(filename,varname,varargin{:});
catch
    % replace missing_value char with a number that can be read
    missing_string=nc_attget(filename,varname,'missing_value');
    Info=nc_getvarinfo(filename,varname);
    if ischar(missing_string) && Info.Nctype~=nc_char
        switch Info.Nctype % cast missing_value into variable's type
            case nc_byte
                missing_value=int8(sscanf(missing_string,'%i'));
            case nc_short
                missing_value=int16(sscanf(missing_string,'%i'));
            case nc_int
                missing_value=int32(sscanf(missing_string,'%i'));
            case nc_float
                missing_value=single(sscanf(missing_string,'%f'));
            case nc_double
                missing_value=double(sscanf(missing_string,'%f'));
        end
    nc_attput(filename,varname,'missing_string',missing_string);
    nc_attput(filename,varname,'missing_value',missing_value);
    end
    result=nc_varget(filename,varname,varargin{:});
end