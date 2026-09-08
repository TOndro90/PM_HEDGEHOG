####################################################################################################
#                             Cosmological Particle Mesh Code HEDGEHOG                             #
####################################################################################################
using HDF5                                                              # HDF5 package for IO
using FFTW                                                              # FFTW package
current_dir = pwd()                                                     # Path to the code
push!(LOAD_PATH, current_dir)
using functions_CIC
using gravity_functions_FFT
using particles_dynamics
using io_functions

############################################# INPUT ################################################
config_file = ARGS[1]                                                   # Get the first command-line argument
include(config_file)                                                    # This will execute the key=value lines in simulation.txt
println("ICs file: ", icsFile)
println("Outputs: ", outputs)
println("Output directory: ", outputDir)
const L_box_const = L_box
const n_Cells_const = n_Cells
const Omega_M_const = Omega_M
const Omega_L_const = Omega_L
const Omega_K_const = Omega_K

####################################################################################################
const Lx = L_box
const Ly = L_box
const Lz = L_box                                
const nx = n_Cells
const ny = n_Cells
const nz = n_Cells
const dx = L_box / n_Cells
const dy = L_box / n_Cells
const dz = L_box / n_Cells
const cosmo_h = H/100
const H0 = H/1000                                                       # [km/s / kpc]
const cosmo_G = 4.30087630e-06                                          # kpc km^2 s^-2 Msun^-1
const xMin = 0
const yMin = 0
const zMin = 0
const xMax = L_box
const yMax = L_box
const zMax = L_box

########################################### LOAD FILES #############################################
#Load ICs
println( "Loading Ics: $(icsFile)")
icsFile = h5open( icsFile, "r")
current_a = read( icsFile, "current_a")
current_z = read( icsFile, "current_z")
mass = read( icsFile, "particle_mass")
pos_x = read( icsFile, "pos_x")
pos_y = read( icsFile, "pos_y")
pos_z = read( icsFile, "pos_z")
vel_x = read( icsFile, "vel_x")
vel_y = read( icsFile, "vel_y")
vel_z = read( icsFile, "vel_z")
const n_local = length(vel_z)
const n_particles = length(vel_z)
close(icsFile)

#################################### MEMORY INITIALIZATION #########################################
density = zeros( Float64, nx, ny, nz )
grav_x = zeros( Float64, n_particles )
grav_y = zeros( Float64, n_particles )
grav_z = zeros( Float64, n_particles )
nGHST_CIC = 1
const nCells_x = nx + 2*nGHST_CIC           
const nCells_y = ny + 2*nGHST_CIC
const nCells_z = nz + 2*nGHST_CIC
gravity_x = zeros( Float64, nz+2*nGHST_CIC, ny+2*nGHST_CIC, nx+2*nGHST_CIC )
gravity_y = zeros( Float64, nz+2*nGHST_CIC, ny+2*nGHST_CIC, nx+2*nGHST_CIC )
gravity_z = zeros( Float64, nz+2*nGHST_CIC, ny+2*nGHST_CIC, nx+2*nGHST_CIC )
nGHST_POT = 3
potential = zeros( Float64, nz+2*nGHST_POT, ny+2*nGHST_POT, nx+2*nGHST_POT)
n_outputs = length( outputs )

########################################### SIMULATION #############################################
println(raw"██╗  ██╗███████╗██████╗  ██████╗ ███████╗██╗  ██╗ ██████╗  ██████╗ ")
println(raw"██║  ██║██╔════╝██╔══██╗██╔════╝ ██╔════╝██║  ██║██╔═══██╗██╔════╝ ")
println(raw"███████║█████╗  ██║  ██║██║  ███╗█████╗  ███████║██║   ██║██║  ███╗")
println(raw"██╔══██║██╔══╝  ██║  ██║██║   ██║██╔══╝  ██╔══██║██║   ██║██║   ██║")
println(raw"██║  ██║███████╗██████╔╝╚██████╔╝███████╗██║  ██║╚██████╔╝╚██████╔╝")
println(raw"╚═╝  ╚═╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ")
println("")
println(raw"This is Cosmological PM code 'HEDGEHOG', version 1.5")
println("")
println("")

println( "Initializing FFT")
G, fft_plan_fwd, fft_plan_bkwd = initialize_FFT( nx, ny, nz, Lx, Ly, Lz, dx, dy, dz )

#Get density
print( "Computing Density CIC")
time_density = @elapsed density =  get_density_CIC( n_local, mass, pos_x, pos_y, pos_z, nCells_x, nCells_y, nCells_z, dx, dy, dz, xMin, yMin, zMin )
dens_avrg = sum( density ) / length( density )
print( "   $(time_density)    $(dens_avrg)\n")

#Get Potential
print( "Computing Potential")
time_potential = @elapsed  get_potential( nx, ny, nz, nGHST_POT, cosmo_G, density, G, fft_plan_fwd, fft_plan_bkwd, potential, dens_avrg, current_a )
print( "   $(time_potential)\n")

#Get Gravitational Field
print( "Computing Gravity Field")
time_grav = @elapsed  get_gravity( nx, ny, nz, nGHST_CIC, nGHST_POT, dx, dy, dz, potential, gravity_x, gravity_y, gravity_z )
print( "   $(time_grav)\n")

#Get gravity CIC accelerations
print( "Computing Gravity CIC")
time_grav_CIC = @elapsed get_gravity_CIC( n_local, pos_x, pos_y, pos_z, xMin, yMin, zMin, dx, dy, dz, gravity_x, gravity_y, gravity_z, grav_x, grav_y, grav_z  )
print( "   $(time_grav_CIC)\n")

#Save Snapshot
n_file = 0
next_output = outputs[1]


println( "\nStarting Simulation")
while current_a < outputs[n_outputs]

  delta_a = get_delta_a_PM(
      n_local,
      current_a,
      vel_x,
      vel_y,
      vel_z,
      dx,
      dy,
      dz,
      Omega_M,
      Omega_L,
      Omega_K,
      H0,
      cosmo_h;
      alpha_DM = 0.3,
      alpha_a = 0.01
  )

  if ( (current_a + delta_a) >= next_output )
    delta_a = (next_output - current_a)
  end

  println( "\nCurrent_a: $(current_a)   delta_a: $(delta_a)"   )


  print( " Updating Particles step 1")
  time_particles_update = @elapsed update_particles_KDK_step_1( 
        n_local, 
        current_a, 
        delta_a, 
        pos_x, 
        pos_y, 
        pos_z, 
        vel_x, 
        vel_y, 
        vel_z, 
        grav_x, 
        grav_y, 
        grav_z, 
        Omega_M, 
        Omega_L, 
        Omega_K, 
        H0, 
        cosmo_h, 
        Lx, 
        Ly, 
        Lz, 
        xMin, 
        yMin, 
        zMin,
        xMax, 
        yMax, 
        zMax 
    )
  print( "   $(time_particles_update)\n")

  #Get density
  print( " Computing Density CIC")
  time_density = @elapsed density =  get_density_CIC( n_local, mass, pos_x, pos_y, pos_z, nCells_x, nCells_y, nCells_z, dx, dy, dz, xMin, yMin, zMin )
  dens_avrg = sum( density ) / length( density )
  print( "   $(time_density)    $(dens_avrg)\n")

  #Get Potential
  print( " Computing Potential")
  time_potential = @elapsed  get_potential( nx, ny, nz, nGHST_POT, cosmo_G, density, G, fft_plan_fwd, fft_plan_bkwd, potential, dens_avrg, current_a+delta_a )
  print( "   $(time_potential)\n")

  #Get Gravitational Field
  print( " Computing Gravity Field")
  time_grav = @elapsed  get_gravity( nx, ny, nz, nGHST_CIC, nGHST_POT, dx, dy, dz, potential, gravity_x, gravity_y, gravity_z )
  print( "   $(time_grav)\n")

  #Get gravity CIC accelerations
  print( " Computing Gravity CIC")
  time_grav_CIC = @elapsed get_gravity_CIC( n_local, pos_x, pos_y, pos_z, xMin, yMin, zMin, dx, dy, dz, gravity_x, gravity_y, gravity_z, grav_x, grav_y, grav_z  )
  print( "   $(time_grav_CIC)\n")

  print( " Updating Particles step 2")
  time_particles_update = @elapsed update_particles_KDK_step_2( n_local, current_a, delta_a, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, grav_x, grav_y, grav_z, Omega_M, Omega_L, Omega_K, H0, cosmo_h , Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax )
  print( "   $(time_particles_update)\n")


  global current_a += delta_a

  if abs(current_a - next_output) < 1e-7 
    global next_output, n_file = write_snapshot(  n_file,  outputDir, current_a, n_local, density,  pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, outputs )
    if abs(current_a - outputs[n_outputs]) < 1e-7
        break
    end
    end
  
end

println( "\nSimulation Finished :)")