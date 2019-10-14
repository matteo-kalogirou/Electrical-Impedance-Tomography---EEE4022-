function    comma2point_overwrite( filespec )
    file    = memmapfile( filespec, 'writable', true );
    comma   = uint8(',');
    point   = uint8('.');
    file.Data( transpose( file.Data==comma) ) = point;
    %%%delete(file)
end