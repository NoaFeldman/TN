@echo off

 matlab -automation -r "prepareEjima; pause(5); quit" > matlab_output.log
 exit
