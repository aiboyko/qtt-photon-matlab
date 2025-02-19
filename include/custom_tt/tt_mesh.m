function mesh=tt_mesh(tts)
tts=flip(tts)
    NT=numel(tts);
    l={};

    for i = 1:NT
        l{i}=tt_ones(2,tts{i}.d);
    end
    
    mesh={};
    for i = 1:NT
        buf=l;
        buf{i}=tts{i};
        mesh{i}=mtkron(buf);
    end
end