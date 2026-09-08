module functions_CIC

export get_density_CIC, get_gravity_CIC
function get_indxs_CIC( n_particles, dx, dy, dz, p_pos_x, p_pos_y, p_pos_z , xMin, yMin, zMin)
  idxs_x = zeros( Int, n_particles )
  idxs_y = zeros( Int, n_particles )
  idxs_z = zeros( Int, n_particles )
  for i in 1:n_particles
      idxs_x[i] = Int(floor( (p_pos_x[i] - 0.5*dx)/dx ))
      idxs_y[i] = Int(floor( (p_pos_y[i] - 0.5*dy)/dy ))
      idxs_z[i] = Int(floor( (p_pos_z[i] - 0.5*dz)/dz ))
  end
  return [idxs_x, idxs_y, idxs_z]
end

function get_density_CIC( n_local, p_mass, p_pos_x, p_pos_y, p_pos_z, nCells_x, nCells_y, nCells_z, dx, dy, dz, xMin, yMin, zMin )
  idxs_x, idxs_y, idxs_z = get_indxs_CIC( n_local, dx, dy, dz, p_pos_x, p_pos_y, p_pos_z, xMin, yMin, zMin )
  #Density array
  rho = zeros( nCells_z, nCells_y, nCells_x )

  for i in 1:n_local
    p_inside = true
    idx_x, idx_y, idx_z = idxs_x[i], idxs_y[i], idxs_z[i]
    if idx_x < -1 || idx_x > nCells_x-3
      p_inside = false
    end
    if idx_y < -1 || idx_y > nCells_y-3
      p_inside = false
    end
    if idx_z < -1 || idx_z > nCells_z-3
      p_inside = false
    end
    if !p_inside
      println( "###### ERROR CIC INDEX OUT OF RANGE ####")
      println( " pInex: $(i)")
      println( " pos_x: $(p_pos_x[i])")
      println( " pos_y: $(p_pos_y[i])")
      println( " pos_z: $(p_pos_z[i])")
      continue
    end
    mass  = p_mass[i]
    pos_x = p_pos_x[i]
    pos_y = p_pos_y[i]
    pos_z = p_pos_z[i]
    cell_center_x = xMin + (idx_x)*dx + 0.5*dx;
    cell_center_y = yMin + (idx_y)*dy + 0.5*dy;
    cell_center_z = zMin + (idx_z)*dz + 0.5*dz;
    delta_x = 1 - (pos_x - cell_center_x)/dx
    delta_y = 1 - (pos_y - cell_center_y)/dy
    delta_z = 1 - (pos_z - cell_center_z)/dz
    idx_x += 2
    idx_y += 2
    idx_z += 2
    # println("$(delta_x)    $(delta_y)    $(delta_z)    ")
    rho[idx_z, idx_y, idx_x]       += mass * delta_x * delta_y * delta_z
    rho[idx_z, idx_y, idx_x+1]     += mass * (1-delta_x) * delta_y * delta_z
    rho[idx_z, idx_y+1, idx_x]     += mass * delta_x * (1-delta_y) * delta_z
    rho[idx_z+1, idx_y, idx_x]     += mass * delta_x * delta_y * (1-delta_z)
    rho[idx_z+1, idx_y+1, idx_x]   += mass * delta_x * (1-delta_y) * (1-delta_z)
    rho[idx_z, idx_y+1, idx_x+1]   += mass * (1-delta_x) * (1-delta_y) * delta_z
    rho[idx_z+1, idx_y, idx_x+1]   += mass * (1-delta_x) * delta_y * (1-delta_z)
    rho[idx_z+1, idx_y+1, idx_x+1] += mass * (1-delta_x) * (1-delta_y) * (1-delta_z)
  end
  rho /= dx*dy*dz
  rho[2, :, :] += rho[end, :, :]
  rho[:, 2, :] += rho[:, end, :]
  rho[:, :, 2] += rho[:, :, end]
  rho[end-1, :, :] += rho[1, :, :]
  rho[:, end-1, :] += rho[:, 1, :]
  rho[:, :, end-1] += rho[:, :, 1]
  return rho[2:end-1, 2:end-1, 2:end-1]
end

function get_gravity_CIC( n_local, p_pos_x, p_pos_y, p_pos_z, xMin, yMin, zMin, dx, dy, dz, gravity_x, gravity_y, gravity_z, p_grav_x, p_grav_y, p_grav_z  )
  idxs_x, idxs_y, idxs_z = get_indxs_CIC(n_local,  dx, dy, dz, p_pos_x, p_pos_y, p_pos_z, xMin, yMin, zMin )
  for i in 1:n_local
    idx_x, idx_y, idx_z = idxs_x[i], idxs_y[i], idxs_z[i]
    pos_x = p_pos_x[i]
    pos_y = p_pos_y[i]
    pos_z = p_pos_z[i]
    cell_center_x = xMin + (idx_x)*dx + 0.5*dx;
    cell_center_y = yMin + (idx_y)*dy + 0.5*dy;
    cell_center_z = zMin + (idx_z)*dz + 0.5*dz;
    delta_x = 1 - (pos_x - cell_center_x)/dx
    delta_y = 1 - (pos_y - cell_center_y)/dy
    delta_z = 1 - (pos_z - cell_center_z)/dz
    idx_x += 2
    idx_y += 2
    idx_z += 2
    g_x, g_y, g_z = get_grav_force_CIC_grid( gravity_x, gravity_y, gravity_z, dx, dy, dz, delta_x, delta_y, delta_z, idx_x, idx_y, idx_z )
    p_grav_x[i] = g_x
    p_grav_y[i] = g_y
    p_grav_z[i] = g_z
  end
end

function get_grav_force_CIC_grid( gravity_x, gravity_y, gravity_z, dx, dy, dz, delta_x, delta_y, delta_z, idx_x, idx_y, idx_z )

  #FOR bottom-left cell0
  g_x_bl = gravity_x[idx_z, idx_y, idx_x]
  g_y_bl = gravity_y[idx_z, idx_y, idx_x]
  g_z_bl = gravity_z[idx_z, idx_y, idx_x]

  #FOR bottom-right
  g_x_br = gravity_x[idx_z, idx_y, idx_x+1]
  g_y_br = gravity_y[idx_z, idx_y, idx_x+1]
  g_z_br = gravity_z[idx_z, idx_y, idx_x+1]

  #FOR bottom-up
  g_x_bu = gravity_x[idx_z, idx_y+1, idx_x]
  g_y_bu = gravity_y[idx_z, idx_y+1, idx_x]
  g_z_bu = gravity_z[idx_z, idx_y+1, idx_x]

  #FOR bottom-right-up
  g_x_bru = gravity_x[idx_z, idx_y+1, idx_x+1]
  g_y_bru = gravity_y[idx_z, idx_y+1, idx_x+1]
  g_z_bru = gravity_z[idx_z, idx_y+1, idx_x+1]

  #FOR top-left
  # X Force component
  g_x_tl = gravity_x[idx_z+1, idx_y, idx_x]
  g_y_tl = gravity_y[idx_z+1, idx_y, idx_x]
  g_z_tl = gravity_z[idx_z+1, idx_y, idx_x]

  #FOR top-right
  g_x_tr = gravity_x[idx_z+1, idx_y, idx_x+1]
  g_y_tr = gravity_y[idx_z+1, idx_y, idx_x+1]
  g_z_tr = gravity_z[idx_z+1, idx_y, idx_x+1]

  #FOR top-up
  g_x_tu = gravity_x[idx_z+1, idx_y+1, idx_x]
  g_y_tu = gravity_y[idx_z+1, idx_y+1, idx_x]
  g_z_tu = gravity_z[idx_z+1, idx_y+1, idx_x]


  #FOR top-right-up
  g_x_tru = gravity_x[idx_z+1, idx_y+1, idx_x+1]
  g_y_tru = gravity_y[idx_z+1, idx_y+1, idx_x+1]
  g_z_tru = gravity_z[idx_z+1, idx_y+1, idx_x+1]


  g_x = g_x_bl*delta_x*delta_y*delta_z         + g_x_br*(1-delta_x)*delta_y*delta_z +
        g_x_bu*delta_x*(1-delta_y)*delta_z     + g_x_bru*(1-delta_x)*(1-delta_y)*delta_z +
        g_x_tl*delta_x*delta_y*(1-delta_z)     + g_x_tr*(1-delta_x)*delta_y*(1-delta_z) +
        g_x_tu*delta_x*(1-delta_y)*(1-delta_z) + g_x_tru*(1-delta_x)*(1-delta_y)*(1-delta_z)

  g_y = g_y_bl*delta_x*delta_y*delta_z         + g_y_br*(1-delta_x)*delta_y*delta_z +
        g_y_bu*delta_x*(1-delta_y)*delta_z     + g_y_bru*(1-delta_x)*(1-delta_y)*delta_z +
        g_y_tl*delta_x*delta_y*(1-delta_z)     + g_y_tr*(1-delta_x)*delta_y*(1-delta_z) +
        g_y_tu*delta_x*(1-delta_y)*(1-delta_z) + g_y_tru*(1-delta_x)*(1-delta_y)*(1-delta_z)

  g_z = g_z_bl*delta_x*delta_y*delta_z         + g_z_br*(1-delta_x)*delta_y*delta_z +
        g_z_bu*delta_x*(1-delta_y)*delta_z     + g_z_bru*(1-delta_x)*(1-delta_y)*delta_z +
        g_z_tl*delta_x*delta_y*(1-delta_z)     + g_z_tr*(1-delta_x)*delta_y*(1-delta_z) +
        g_z_tu*delta_x*(1-delta_y)*(1-delta_z) + g_z_tru*(1-delta_x)*(1-delta_y)*(1-delta_z)

  return [ g_x g_y g_z ]
end


end