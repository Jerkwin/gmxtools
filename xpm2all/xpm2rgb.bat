@ echo off

set File=%~n1

bash xpm2all.bsh -rgb 8 %File%

start /min D:\PicTool\iview-4.54\i_view64.exe %File%~iv.xpm

sleep 1
del /F /Q %File%~iv.xpm
