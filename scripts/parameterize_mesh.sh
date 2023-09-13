#! /bin/bash

bin_dir=../build/bin
input_dir=$1
fname=$2
output_dir=$3

${bin_dir}/generate_hidden_cones -i ${input_dir} -f ${fname} -o ${output_dir}
# Note: Cone refinement requires python libraries, including libigl
#python3 refine_cone_neighborhoods.py -i ${output_dir} -f ${fname} -o ${output_dir}
${bin_dir}/parameterize_mesh -i ${output_dir} -f ${fname} -o ${output_dir}
