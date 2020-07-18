ch1=classHelloHandle(1)
ch2=classHelloHandle(2)
a=[ch1,ch2]
a(1).veshalka=11
a(2).veshalka=22
ch1
ch2

clear ch1
a(2)
%мораль: clear удаляет именно данную ссылку, а не ссылку и объект в памяти
% destroy('ch1')
% a(2)