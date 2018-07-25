@echo off
X:
cd zzhang\Python\pyRTAOI
call activate pyRTAOI
call pyuic5 GUI.ui -o GUI.py
echo converted 
pause