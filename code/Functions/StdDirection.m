function v_new = StdDirection(v_old)

v_new = v_old;

[~, max_dir] = max(abs(v_old)); %1 -> x, 2 -> y, 3 -> z

if(v_old(max_dir)<0)
    
    v_new = -v_old;
    
end

end