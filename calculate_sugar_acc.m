function [efficiency,growth_capacity]=calculate_sugar_acc(t,p,eme,efficiency,growth_capacity,is_light)
if t < eme + 24 + 1
    efficiency(t) = p(84);
    growth_capacity(t) = p(63);
else

    growth_capacity(t) = p(63);


    if is_light(t) == 1 && is_light(t-1) == 0

        efficiency(t) = p(84);
    else
        efficiency(t) = efficiency(t-1);
    end
end    