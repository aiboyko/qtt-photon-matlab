function y= oldcell2core(ttv_oldcells)
p=@(x)permute(x,[2 1 3]);
ttv_newcells=cellfun(p,ttv_oldcells,'UniformOutput',false);
y=cell2core(tt_tensor,ttv_newcells);
end