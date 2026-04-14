@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="bklab-27" (taskkill /f /pid 14604)
if /i "%LOCALHOST%"=="bklab-27" (taskkill /f /pid 27396)
if /i "%LOCALHOST%"=="bklab-27" (taskkill /f /pid 29920)
if /i "%LOCALHOST%"=="bklab-27" (taskkill /f /pid 6192)

del /F cleanup-ansys-bklab-27-6192.bat
