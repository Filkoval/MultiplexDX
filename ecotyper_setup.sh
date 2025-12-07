#!/bin/bash
#ECOTYPER SETUP
#This script creates configuration and input file for Ecotyper
#Arguments that are supported by this script:
	# -s  =name of the sample
	# -i  =path to input directory(SpaceRanger output)
	# -o  =(optional)path to output directory, default: sample_name_ecotyper_output
	# -c  =(optional)path to configuration .yml file for ecotyper
	# -r  =(optional)remove temporary files(ecotyper input, config file)
	# -e  =(optional)run ecotyper on the created files

remove_temps=false
ecotyper=false
while getopts "s:i:o:c:re" opt; do
    case $opt in
	s) sample_name="$OPTARG";;
        i) input_dir="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
	c) config_file="$OPTARG" ;;
	r) remove_temps=true ;;
	e) ecotyper=true ;;
        *) echo "Incorrect flag use, supported flags are: -s, -i, -o, -c, -r"; exit 1 ;;
    esac
done


if [[ -d "$output_dir" ]]; then
    rm -r "$output_dir"
fi

if [[ ! -f EcoTyper_recovery_visium.R ]]; then
    echo "Please move this script into ecotyper repo"
    exit 1
fi

if [[ -z "$sample_name" || -z "$input_dir" ]]; then
	echo 'Please add all mandatory flags with correct argument: -s "sample_name" -i "input_directory(SpaceRanger output directory)"'
	exit 1
fi

if [[ -n "$sample_name" && -z "$output_dir" ]]; then
    output_dir="${sample_name}_ecotyper_output"
fi

if [[ -n "$sample_name" && -z "$config_file" ]]; then
    config_file="${sample_name}_config.yml"
fi

extension="${config_file##*.}"

if [[ "$extension" != "yml" ]]; then
    echo "config file does not have a correct format - format should be filename.yml"
fi

mkdir -p $output_dir/temp

echo "Output directory $output_dir was created"

discovery_dataset="Carcinoma"
recovery_dataset="VisiumBreast"
recovery_cell_types="NULL"
background="Epithelial.cells"
username="filkoval@natur.cuni.cz"
token="d104a21f780b61cf1c6831c1e5f91e80"
threads=10
sing_path=NULL

#Create ecotyper input directory

cp $input_dir/filtered_feature_bc_matrix/* $output_dir/temp/
cp $input_dir/spatial/tissue_positions.csv $output_dir/temp/tissue_positions_list.csv


#Create config file


cat <<EOF > "$config_file"
default:
  Input:
    Discovery dataset name: "$discovery_dataset"
    Recovery dataset name: "$recovery_dataset"
    Input Visium directory: "$output_dir/temp"
    Recovery cell type fractions: "$recovery_cell_types"
    Background cell type: "$background"
    CIBERSORTx username: "$username"
    CIBERSORTx token: "$token"

  Output:
    Output folder: "$output_dir"

  Pipeline settings:
    Number of threads: $threads
    CIBERSORTx fractions Singularity path: $sing_path
EOF

echo "Configuration file $config_file was created"

if [[ "$ecotyper" == "true" ]]; then
	Rscript EcoTyper_recovery_visium.R -c "$config_file"
fi


if [[ "$remove_temps" == "true" ]]; then
	rm -r $output_dir/temp/
	rm $config_file
fi
