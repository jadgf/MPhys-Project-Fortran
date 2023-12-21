input_file_path = '2d_trajectory.dat'
output_file_path = 'clean_2d_trajectory.dat'

# Open the input file in read mode
with open(input_file_path, 'r') as input_file:
    # Read all lines from the file
    lines = input_file.readlines()

# Filter out lines with only zeroes
non_zero_lines = [line for line in lines if any(float(value) != 0 for value in line.split())]

# Open the output file in write mode
with open(output_file_path, 'w') as output_file:
    # Write the modified lines back to the file
    output_file.writelines(non_zero_lines)

print(f"Lines with only zeroes removed. Modified content saved to {output_file_path}")

