@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="bklab-25" (taskkill /f /pid 16192)
if /i "%LOCALHOST%"=="bklab-25" (taskkill /f /pid 17592)
if /i "%LOCALHOST%"=="bklab-25" (taskkill /f /pid 6760)
if /i "%LOCALHOST%"=="bklab-25" (taskkill /f /pid 23368)

del /F cleanup-ansys-bklab-25-23368.bat
