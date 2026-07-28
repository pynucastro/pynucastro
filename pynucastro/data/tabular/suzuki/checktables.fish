for f in (ls *-toki);
    if test "FAILED" = (./checktables.py $f);
        echo failed;
    end
end
