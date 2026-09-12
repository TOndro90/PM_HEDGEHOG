module particles_dynamics

export  update_particles_KDK_step_1,
        update_particles_KDK_step_2,
        get_delta_a_PM

########################################
function adot( 
    a, 
    Omega_M, 
    Omega_L, 
    Omega_K, 
    H0, 
    cosmo_h 
)
  return  H0* a * sqrt( Omega_M/a^3 + Omega_L + Omega_K/a^2 ) / cosmo_h 
end

########################################
function get_delta_a_PM(
    n_local,
    current_a,
    p_vel_x,
    p_vel_y,
    p_vel_z,
    dx,
    dy,
    dz,
    omega_M,
    omega_L,
    omega_K,
    H0,
    cosmo_h;
    alpha_DM, 
    alpha_a
)
    vmax_x = 0.0
    vmax_y = 0.0
    vmax_z = 0.0
    for i in 1:n_local
        vx = abs(p_vel_x[i])
        vy = abs(p_vel_y[i])
        vz = abs(p_vel_z[i])
        if vx > vmax_x
            vmax_x = vx
        end
        if vy > vmax_y
            vmax_y = vy
        end
        if vz > vmax_z
            vmax_z = vz
        end
    end
    if vmax_x > 0.0
        delta_t_x = alpha_DM * current_a * dx / vmax_x
    else
        delta_t_x = Inf
    end
    if vmax_y > 0.0
        delta_t_y = alpha_DM * current_a * dy / vmax_y
    else
        delta_t_y = Inf
    end
    if vmax_z > 0.0
        delta_t_z = alpha_DM * current_a * dz / vmax_z
    else
        delta_t_z = Inf
    end
    delta_t = min( delta_t_x, delta_t_y, delta_t_z )
    adot_now = adot( current_a, omega_M, omega_L, omega_K, H0, cosmo_h )
    delta_a_DM = adot_now * delta_t
    delta_a_expansion = alpha_a * current_a
    delta_a = min( delta_a_DM, delta_a_expansion )
    return delta_a
end

########################################
function update_particles_KDK_step_1( n_local, current_a, delta_a, p_pos_x, p_pos_y, p_pos_z, p_vel_x, p_vel_y, p_vel_z, p_grav_x,  p_grav_y, p_grav_z, omega_M, omega_L, omega_K, H0, cosmo_h, Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax  )

  a_half = current_a + delta_a/2
  delta_t = delta_a / adot( current_a, omega_M, omega_L, omega_K, H0, cosmo_h )

  for i in 1:n_local
    p_vel_x[i] = p_vel_x[i] * current_a / a_half +  p_grav_x[i] *delta_t / ( 2*a_half ) 
    p_vel_y[i] = p_vel_y[i] * current_a / a_half +  p_grav_y[i] *delta_t / ( 2*a_half ) 
    p_vel_z[i] = p_vel_z[i] * current_a / a_half +  p_grav_z[i] *delta_t / ( 2*a_half ) 

    p_pos_x[i] += p_vel_x[i] * delta_t / a_half 
    p_pos_y[i] += p_vel_y[i] * delta_t / a_half 
    p_pos_z[i] += p_vel_z[i] * delta_t / a_half 
    
    if p_pos_x[i] >= xMax
      p_pos_x[i] -= Lx
    end
    if p_pos_x[i] < xMin
      p_pos_x[i] += Lx
    end

    if p_pos_y[i] >= yMax
      p_pos_y[i] -= Ly
    end
    if p_pos_y[i] < yMin
      p_pos_y[i] += Ly
    end

    if p_pos_z[i] >= zMax
      p_pos_z[i] -= Lz
    end
    if p_pos_z[i] < zMin
      p_pos_z[i] += Lz
    end

  end

end

########################################
function update_particles_KDK_step_2( n_local, current_a, delta_a, p_pos_x, p_pos_y, p_pos_z, p_vel_x, p_vel_y, p_vel_z, p_grav_x,  p_grav_y, p_grav_z, omega_M, omega_L, omega_K, H0, cosmo_h, Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax  )
  a_half = current_a + delta_a/2
  an = current_a + delta_a
  delta_t = delta_a / adot( current_a, omega_M, omega_L, omega_K, H0, cosmo_h )
  for i in 1:n_local
    p_vel_x[i] = p_vel_x[i] * a_half / an +  p_grav_x[i] *delta_t / ( 2*an ) 
    p_vel_y[i] = p_vel_y[i] * a_half / an +  p_grav_y[i] *delta_t / ( 2*an ) 
    p_vel_z[i] = p_vel_z[i] * a_half / an +  p_grav_z[i] *delta_t / ( 2*an ) 
  end

end


end