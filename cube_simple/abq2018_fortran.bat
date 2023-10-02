@echo off
setlocal
set ABA_COMMAND=%~nx0
set ABA_COMMAND_FULL=%~f0
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat" intel64 vs2019
"C:\SIMULIA\CAE\2018\win_b64\code\bin\ABQLauncher.exe" %*
endlocal