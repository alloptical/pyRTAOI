function stars = get_stars_for_P (P)
if P<0.0001
    stars = '****';
else if P<0.001
        stars = '***';
    else if(P<0.01)
            stars = '**';
        else if (P<0.05)
                stars = '*';
            else
                stars = 'n.s.';
            end
        end
    end
end
end

