#! /bin/bash

output_dir=output
bin_dir=../build/bin
data_dir=../data/meshes
camera_dir=../data/cameras

model_list=(\
 "bigguy-cleaned-2"\
 "blub_quadrangulated_tri_clean"\
 "bob_quad_tri_clean"\
 "fertility_tri_clean"\
 "killeroo_small_tri_clean"\
 "monsterfrog-cleaned-2"\
 "ogre_controlmesh-cleaned-3"\
 "pawn_tri_clean"\
 "pipes_tri_clean"\
 "spot_quadrangulated_tri_clean"\
)

for model in ${model_list[@]}
do
  mkdir -p ${output_dir}/fig-example/${model}
  ${bin_dir}/generate_example_figure \
    -i ${data_dir}/${model}_conf_simplified_with_uv.obj \
    -c ${camera_dir}/${model}_camera.txt \
    -o ${output_dir}/fig-example/${model}
done
