t=5;
tic;
for i=1:1e6;
a=mod(5,3);
end
toc

tic
f=@(x) x>4;

for i=1:1e6;
    if a>3
        a=1;
    end
end
toc

tic
for i=1:1e6;

        a=1;

end
toc