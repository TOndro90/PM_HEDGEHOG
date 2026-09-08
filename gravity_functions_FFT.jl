module gravity_functions_FFT
using FFTW

export get_potential, initialize_FFT, get_gravity
# Compute the potential from the input density

function get_potential( nx, ny, nz, nGHST, Gconst, density, G, fft_plan_f, fft_plan_b, potential, dens_avrg, current_a )
    rho_trans = fft_plan_f * (4*pi*Gconst * (density .- dens_avrg) ) /current_a
    rho_trans = G .* rho_trans
    rho_trans[1] = 0
    # Apply inverse FFTW
    phi = real( fft_plan_b * rho_trans )
    potential[nGHST+1:end-nGHST, nGHST+1:end-nGHST, nGHST+1:end-nGHST] = phi
    potential[1, :, : ] = potential[end-2*nGHST+1,:,:]
    potential[2, :, : ] = potential[end-2*nGHST+2,:,:]
    potential[3, :, : ] = potential[end-2*nGHST+3,:,:]

    potential[:, 1, : ] = potential[:,end-2*nGHST+1,:]
    potential[:, 2, : ] = potential[:,end-2*nGHST+2,:]
    potential[:, 3, : ] = potential[:,end-2*nGHST+3,:]
    
    potential[:, :, 1 ] = potential[:,:,end-2*nGHST+1]
    potential[:, :, 2 ] = potential[:,:,end-2*nGHST+2]
    potential[:, :, 3 ] = potential[:,:,end-2*nGHST+3]

    potential[end,:,:] = potential[nGHST+3,:,:]
    potential[end-1,:,:] = potential[nGHST+2,:,:]
    potential[end-2,:,:] = potential[nGHST+1,:,:]

    potential[:,end,:] = potential[:,nGHST+3,:]
    potential[:,end-1,:] = potential[:,nGHST+2,:]
    potential[:,end-2,:] = potential[:,nGHST+1,:]    
    
    potential[:,:,end] = potential[:,:,nGHST+3]
    potential[:,:,end-1] = potential[:,:,nGHST+2]
    potential[:,:,end-2] = potential[:,:,nGHST+1] 
    return potential
  end

function initialize_FFT( nCells_x, nCells_y, nCells_z, Lx, Ly, Lz, dx, dy, dz )
  dens = zeros( Float64,nCells_z, nCells_y, nCells_x )
  # Initialize FFTW
  fft_plan_fwd  = plan_fft( dens, flags=FFTW.ESTIMATE, timelimit=20 )
  fft_plan_bkwd = plan_ifft( dens, flags=FFTW.ESTIMATE, timelimit=20 )

  fft_kx = sin.( pi/Lx*dx * range(0, nCells_x-1, nCells_x) ).^2
  fft_ky = sin.( pi/Ly*dy * range(0, nCells_y-1, nCells_y) ).^2
  fft_kz = sin.( pi/Lz*dz * range(0, nCells_z-1, nCells_z) ).^2
  G = zeros(Float64, nCells_z, nCells_y, nCells_x )
  for i in 1:nCells_x
    sin_kx = fft_kx[i]
    for j in 1:nCells_y
      sin_ky = fft_ky[j]
      for k in 1:nCells_z
        sin_kz =  fft_kz[k]
        G[k,j,i] = -1/ ( sin_kx + sin_ky + sin_kz )  *dx^2/4  #NOTE: DX^2 only for dx=dy=dz
      end
    end
  end
  G[1,1,1] = 0
  return [ G, fft_plan_fwd, fft_plan_bkwd ]
end

        
  function get_gravity( nx, ny, nz, nGHST_CIC, nGHST_POT, dx, dy, dz, potential, grav_x, grav_y, grav_z )
    nx_g = nx+2*nGHST_CIC
    ny_g = ny+2*nGHST_CIC
    nz_g = nz+2*nGHST_CIC
    for i in 1:nx_g
      for j in 1:ny_g
        for k in 1:nz_g
          i_p, j_p, k_p = i+2, j+2, k+2
          pot2_l = potential[ k_p, j_p, i_p-2 ]
          pot1_l = potential[ k_p, j_p, i_p-1 ]
          pot1_r = potential[ k_p, j_p, i_p+1 ]
          pot2_r = potential[ k_p, j_p, i_p+2 ]
          grav_x[k, j, i] = - (pot2_l - 8*pot1_l + 8*pot1_r - pot2_r ) / (12*dx)
        end
      end
    end
  
    for i in 1:nx_g
      for j in 1:ny_g
        for k in 1:nz_g
          i_p, j_p, k_p = i+2, j+2, k+2
          pot2_l = potential[ k_p, j_p-2, i_p ]
          pot1_l = potential[ k_p, j_p-1, i_p ]
          pot1_r = potential[ k_p, j_p+1, i_p ]
          pot2_r = potential[ k_p, j_p+2, i_p ]
          grav_y[k, j, i] = - (pot2_l - 8*pot1_l + 8*pot1_r - pot2_r ) / (12*dy)
        end
      end
    end
  
    for i in 1:nx_g
      for j in 1:ny_g
        for k in 1:nz_g
          i_p, j_p, k_p = i+2, j+2, k+2
          pot2_l = potential[ k_p-2, j_p, i_p ]
          pot1_l = potential[ k_p-1, j_p, i_p ]
          pot1_r = potential[ k_p+1, j_p, i_p ]
          pot2_r = potential[ k_p+2, j_p, i_p ]
          grav_z[k, j, i] = - (pot2_l - 8*pot1_l + 8*pot1_r - pot2_r ) / (12*dz)
        end
      end
    end
    return grav_x, grav_y, grav_z
  end

end