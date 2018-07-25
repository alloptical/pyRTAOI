@echo off
set root=C:\Users\Zihui\Anaconda3
call %root%\Scripts\activate.bat %root%
X:
cd Python\pyRTAOI
call activate C:\Users\Zihui\Anaconda3\envs\pyRTAOI
call pyuic5 GUI.ui -o GUI.py
echo converted 
pause