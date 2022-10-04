export initialize

function initialize(configurations::ETOInitialDataConfigurations)

    
    packets = my_zeros(configurations)

    index = 1
    for initial_position in get_initial_positions(configurations)

    for index in range(number_of_points)

            @views packets_at_position = packets[:, index]
            initialize_single!(ray, initial_time, pixel_coordinates, configurations, container)
    end

    return rays

end

function initialize_packets_at_position!(packets, position, configurations)

    number_of_packets = size(positions,2)

    N = N_ep * N_vpp   

    vn = np.zeros((N,4), dtype=dtype)
    pe = np.zeros((N,4), dtype=dtype)

    for j in range(N_ep):
        v = generate_halfsphere_vector(N_vpp)
        gj = g(p[j,:],M,A)
        e = generate_local_frames(p[j,:],t[j,:],gj)
        for i in range(N_vpp):   
            
            v[i,:] = e[0,:] + e[1,:] * v[i,1]  + e[2,:] * v[i,2] + e[3,:] * v[i,3]
            vn[N_vpp * j + i,:] = v[i,:]
            pe[N_vpp * j + i,:] = p[j,:]
    
        end
    end
   
    return np.block([pe, vn])
end

