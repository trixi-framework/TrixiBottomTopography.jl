
# Converting DGM data into readable files (1D)
"""
    convert_dgm_1d(path_read::String, path_write::String; excerpt = 1, direction = "x", section = 1)

Function to convert [DGM](https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/)
data files into one dimensional readable files.

Inputs:
- `path_read`: String of the path of the DGM data which should be converted
- `path_write`: String of the path where the new file shuold be saved.
              (Needs to also include the name of the file)
- `excerpt`: Optional integer that specifies a stride through of the data that will be extracted. E.g.
           if excerpt is set to 10, only every 10th `x` and `y` value are considered with their
           corresponsing `z` values. The default value is 1 which means that every value is taken.
- `direction`: Optional String that specifies if the one dimensional data should be read from the
             `x` or `y` direction. By default this is set to the `x` direction.
- `section`: Optional interger which can be between 1 and 1000 and specifies which section of the
           other dimension should be chosen. By default this values is set to 1 which means
           that if direction is set to `x`, the corresponding `z` values are taken with respect to the
           first `y` value.
"""
function convert_dgm_1d(path_read::String, path_write::String;
                        excerpt = 1, direction = "x", section = 1)

  # Check if section has an accepted value
  if (section < 1 | section > 1000)
    @error("The value for section must be between 1 and 1000!")
  end

  # Get file data
  read_file = open(path_read)
  data = readlines(read_file)
  close(read_file)

  # Get dimensions
  length_data    = length(data)
  dimension_data = Int(sqrt(length_data))

  if direction == "x"

    # Create interim vectors which save all x and z values
    # accordingly
    x_all = zeros(length_data)
    z_all = zeros(length_data)

    # Save the values in the corresponding vectors
    for i in 1:length_data

      line = split(data[i], " ")

      x_all[i] = parse(Float64, line[1])
      z_all[i] = parse(Float64, line[3])

    end

    # Save only the relevant data
    x_uniq = x_all[1:excerpt:dimension_data]
    y_uniq = reshape(z_all, (dimension_data, dimension_data))[1:excerpt:dimension_data, section]

  elseif direction == "y"

    # Create interim vectors which save all x and z values
    # accordingly
    y_all = zeros(length_data)
    z_all = zeros(length_data)

    # Save the values in the corresponding vectors
    for i in 1:length_data

      line = split(data[i], " ")

      y_all[i] = parse(Float64, line[2])
      z_all[i] = parse(Float64, line[3])

    end

    # Save only the relevant data
    x_uniq = y_all[1:(excerpt*dimension_data):length_data]
    y_uniq = reshape(z_all, (dimension_data, dimension_data))[section, 1:excerpt:dimension_data]

  else
    @error("The input direction can either be the String x or y")
  end

  # Write the data to the file
  write_file = open(path_write, "w")

  # The resulting file has the following form:
  #
  # # Number of x values
  # Number of x values
  # # x values
  # x_1
  # ...
  # x_n
  # # y values
  # y_1
  # ...
  # y_n

  println(write_file, "# Number of x values")
  println(write_file, "$(Int(floor(dimension_data/excerpt)))")
  println(write_file, "# x values")
  for i in eachindex(x_uniq)
    println(write_file, x_uniq[i])
  end
  println(write_file, "# y values")
  for i in eachindex(y_uniq)
    println(write_file, y_uniq[i])
  end

 close(write_file)
end


# Converting DGM data into readable files (2D)
"""
    convert_dgm_2d(path_read::String, path_write::String; excerpt = 1)

Function to convert [DGM](https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/)
data files into two dimensional readable files.

Inputs:
- `path_read`: String of the path of the DGM data which should be converted
- `path_write`: String of the path where the new file shuold be saved.
              (Needs to also include the name of the file)
- `excerpt`: Optional integer that specifies a stride through of the data that will be extracted. E.g.
           if excerpt is set to 10, only every 10th `x` and `y` value are considered with their
           corresponsing `z` values. The default value is 1 which means that every value is taken.
"""
function convert_dgm_2d(path_read::String, path_write::String; excerpt = 1)

  # Get file data
  read_file = open(path_read)
  data = readlines(read_file)
  close(read_file)

  # Get dimensions
  length_data    = length(data)
  dimension_data = Int(sqrt(length_data))

  # Create interim vectors which save all x, y and z values accordingly
  x_all = zeros(length_data)
  y_all = zeros(length_data)
  z_all = zeros(length_data)

  # Save the values in the corresponding vectors
  for i in 1:length_data

    line = split(data[i], " ")

    x_all[i] = parse(Float64, line[1])
    y_all[i] = parse(Float64, line[2])
    z_all[i] = parse(Float64, line[3])

  end

  # Save only the relevant data
  x_uniq = x_all[1:excerpt:dimension_data]
  y_uniq = y_all[1:(excerpt*dimension_data):length_data]
  z_uniq = reshape(z_all, (dimension_data, dimension_data)
                  )[1:excerpt:dimension_data,
                    1:excerpt:dimension_data]

  # Write the data to the file
  write_file = open(path_write, "w")

  # The resulting file has the following form:
  #
  # # Number of x values
  # Number of x values
  # # Number of y values
  # Number of y values
  # # x values
  # x_1
  # ...
  # x_n
  # # y values
  # y_1
  # ...
  # y_m
  # # z values
  # z_11
  # z_12
  # ...
  # z_mn

  println(write_file, "# Number of x values")
  println(write_file, "$(Int(floor(dimension_data/excerpt)))")
  println(write_file, "# Number of y values")
  println(write_file, "$(Int(floor(dimension_data/excerpt)))")
  println(write_file, "# x values")
  for i in eachindex(x_uniq)
    println(write_file, x_uniq[i])
  end
  println(write_file, "# y values")
  for i in eachindex(y_uniq)
    println(write_file, y_uniq[i])
  end
  println(write_file, "# z values")
  for i in eachindex(z_uniq)
    println(write_file, z_uniq[i])
  end

  close(write_file)
end
