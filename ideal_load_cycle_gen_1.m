function [ideal_load_cycle_data_1]=ideal_load_cycle_gen_1(strain_rate)

sec_n=100;

hold_check=[0.1,0.2,0.3,0.4,0.5];

turn_point_times=[0.1/strain_rate;
                  9000;
                  0.1/strain_rate;
                  9000;
                  0.1/strain_rate;
                  9000;
                  0.1/strain_rate;
                  9000;
                  0.1/strain_rate;
                  9000;
                  0.5/strain_rate;
                  0.5/strain_rate;
                  0.5/strain_rate];
              
turn_point_times=cumsum(turn_point_times);         

turn_point_strain=[0.1;
                   0.1;
                   0.2;
                   0.2;
                   0.3;
                   0.3;
                   0.4;
                   0.4;
                   0.5;
                   0.5;
                   0;
                   -0.5;
                   0];
               
               
 for ii=1:1:size(turn_point_times)
    if ii==1
        t0=0;
        e0=0;        
        
        ideal_load_cycle_data_1(1:(sec_n+1),1)=transpose(linspace(t0,turn_point_times(ii,1),sec_n+1));
        ideal_load_cycle_data_1(1:(sec_n+1),2)=transpose(linspace(e0,turn_point_strain(ii,1),sec_n+1));
        clear t0 e0
    else
        if turn_point_strain(ii,1)==turn_point_strain(ii-1,1)
            t0=turn_point_times(ii-1,1);
            e0=turn_point_strain(ii-1,1);

            %t_temp=transpose(logspace(log10(t0+1),log10(turn_point_times(ii,1)),sec_n));
            %e_temp=transpose(logspace(log10(e0),log10(turn_point_strain(ii,1)),sec_n));
            t_temp=transpose(linspace(t0,turn_point_times(ii,1),sec_n+1));
            e_temp=transpose(linspace(e0,turn_point_strain(ii,1),sec_n+1));

            ideal_load_cycle_data_1((((ii-1)*sec_n)+2):((ii*sec_n)+1),1)=t_temp(2:end);
            ideal_load_cycle_data_1((((ii-1)*sec_n)+2):((ii*sec_n)+1),2)=e_temp(2:end);
            clear t0 e0 t_temp e_temp
        else        
            t0=turn_point_times(ii-1,1);
            e0=turn_point_strain(ii-1,1);

            t_temp=transpose(linspace(t0,turn_point_times(ii,1),sec_n+1));
            e_temp=transpose(linspace(e0,turn_point_strain(ii,1),sec_n+1));

            ideal_load_cycle_data_1((((ii-1)*sec_n)+2):((ii*sec_n)+1),1)=t_temp(2:end);
            ideal_load_cycle_data_1((((ii-1)*sec_n)+2):((ii*sec_n)+1),2)=e_temp(2:end);
            clear t0 e0 t_temp e_temp
        end
    end
end
clear ii