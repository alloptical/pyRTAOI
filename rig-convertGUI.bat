@echo off
cd C:\Users\User\Desktop\pyRTAOI-rig
call activate caiman
call pyuic5 GUI.ui -o GUI.py
echo converted 
pause