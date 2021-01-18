call "C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\Professional\\VC\\Auxiliary\\Build\\vcvarsall.bat" amd64

set MRH_BASE=%cd%\\..\\..
set GEOGRAPHICLIB_BASE=%MRH_BASE%\\..\\GeographicLib-1.51

set INCLUDE=%INCLUDE%;%MRH_BASE%\\include;%GEOGRAPHICLIB_BASE%\\include;%cd%

set LIB=%LIB%;%MRH_BASE%\\lib;%cd%

set PATH=%PATH%;%MRH_BASE%\\bin;%cd%

"C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\Professional\\Common7\\IDE\\devenv.exe" /USEENV MultiResHealpix.sln