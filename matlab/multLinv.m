    function[X] = multLinv(lfile, lenl, X)
        for i=1:lenl
            ri = lfile(1,i);
            ci = lfile(2,i);
            li = lfile(3,i);

            X(ri,:) =X(ri,:) + li*X(ci,:); 
        end
    end
