    function [X] = multL(lfile, lenl,X)
        for i = lenl:-1:1
            ri = lfile(1,i);
            ci = lfile(2,i);
            li = lfile(3,i);

            X(ri,:) = X(ri,:) - li*X(ci,:);
        end
    end

