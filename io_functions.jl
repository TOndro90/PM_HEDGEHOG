module io_functions

using HDF5

export write_snapshot

function write_snapshot( n_file,  outputDir, current_a, n_local,  density, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, outputs )

  println( "\nWriting Snapshot: $(n_file)    current_a: $(current_a)" )
  fileName = "$(n_file)_particles.h5"
  outFile = h5open( outputDir * fileName, "w")

  current_z = 1/current_a - 1
  attrs(outFile)["current_a"] = current_a
  attrs(outFile)["current_z"] = current_z
  attrs(outFile)["n_particles_local"] = n_local

  outFile[ "density" ] = density

  outFile[ "pos_x" ] = pos_x
  outFile[ "pos_y" ] = pos_y
  outFile[ "pos_z" ] = pos_z

  outFile[ "vel_x" ] = vel_x
  outFile[ "vel_y" ] = vel_y
  outFile[ "vel_z" ] = vel_z

  close(outFile)
  n_file += 1
  println( " Snapshot saved:  $(fileName)")
  if n_file+1 > length( outputs )
    return 0, 0
  end

  next_output = outputs[n_file+1]
  println( " Next Output:  $(next_output)")
  return next_output, n_file
end

end