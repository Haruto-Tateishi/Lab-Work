# this script runs the script SF_Ratios.py automatically for all sfs files that are ready to be analyzed.

import subprocess
import os
import chardet

envir = os.environ.copy()
envir['LC_ALL'] = 'en_US.UTF-8'
envir['LANG'] = 'en_US.UTF-8'

# Path to the script you want to run (replace with your actual script)
script_path = "/home/haruto/Lab-Work/SF_Ratios/SF_Ratios.py"

# Directory containing the input files
input_dir = "/home/haruto/Lab-Work/Document/SFS/FlankingBases/syn_codon/single/1kG"
output_dir = "/home/haruto/Lab-Work/Document/SF_Ratios/FlankingBases/syn_codon/1kG"


# Get a sorted list of all files in the input directory
input_files = sorted([f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))])

i = 0
# Loop through each input file
for idx, input_file in enumerate(input_files):
    input_file_path = os.path.join(input_dir, input_file)

    if i == 0 or i == 1:
      label = "1kG_" + input_file.split(sep="_")[7]
    else:
      label = "1kG_" + input_file.split(sep="_")[6]
    print(label)
    i +=1

    # Construct the command with options -i and -o
    command = ["python3", script_path, "-d", "fixed2Ns", "-f", "unfolded", "-i", "3", "-u", "-p", label, "-a", input_file_path, "-r", output_dir]

    try:
        # Execute the command
        subprocess.run(command, check=True, env=envir)
        print(f"Successfully processed {input_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {input_file}: {e}")

