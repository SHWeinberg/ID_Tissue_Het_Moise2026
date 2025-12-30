function ts = get_time_variable(trange, dt1, dt2,dt1_samp, dt2_samp, twin, bcl)

    ti = trange(1);  
    count = 1;
    ts = [];
    while ti < trange(2)
        if mod(ti, bcl)<twin
            dt = dt1; dt_samp = dt1_samp; 
        else
            dt = dt2; dt_samp = dt2_samp; 
        end
        
        %ts def identical as in main loop
        if ~mod(ti, dt_samp) 
            ts(count) = ti;
            count = count + 1;
        end
    
        % ti = round(ti + dt_samp,5);
        ti = round(ti + dt,5);
    end
    
end