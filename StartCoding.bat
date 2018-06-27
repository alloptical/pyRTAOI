@echo off
set root=C:\Users\Zihui\Anaconda3
call %root%\Scripts\activate.bat %root%
Y:
cd Y:\zzhang\Python\pyRTAOI
call activate C:\Users\Zihui\Anaconda3\envs\pyRTAOI
start spyder
echo ALL SET
exit