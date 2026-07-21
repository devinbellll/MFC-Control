function blkStruct = slblocks
%SLBLOCKS Register the MFC block library with the Simulink Library Browser.
%
%   Simulink calls this automatically for every folder on the MATLAB path.
%   It only takes effect once library/mfc_lib.mdl exists -- generate it with
%   build_mfc_lib -- and after a  >> sl_refresh_customizations.
%
%   See also BUILD_MFC_LIB.

    Browser.Library = 'mfc_lib';
    Browser.Name    = 'Model-Free Control (MFC)';
    blkStruct.Browser = Browser;
end
