@ echo off

set File=%~n1

bash xpm2all.bsh -rgb 8 %File%

start /min C:\Users\Jicun\D_Sft\PicTool\IrfanView\i_view32.exe %File%~iv.xpm

sleep 1
del /F /Q %File%~iv.xpm
