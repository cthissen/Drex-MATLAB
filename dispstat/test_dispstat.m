dispstat('','init'); % One time only initialization 
dispstat(sprintf('\nBegining the process...'),'keepthis','timestamp'); 
for i = 1:100 
dispstat(sprintf('Progress %d%%',i),'timestamp'); 
pause(0.1)
%doing some heavy stuff here 
end 
dispstat('Finished.','keepprev');