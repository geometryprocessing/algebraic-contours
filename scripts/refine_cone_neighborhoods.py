# Script to refine neighborhoods of cone angles with loop subdivision.

import os, argparse
import numpy as np
import igl
import pymeshlab
ms = pymeshlab.MeshSet()

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Script to refine neighborhoods of cones")
    parser.add_argument("-f", "--fname",         help="filename of the obj file")
    parser.add_argument("-i", "--input_dir",     help="input folder that stores obj files and Th_hat")
    parser.add_argument("-o", "--output_dir",     help="output folder for refined mesh")
    args = parser.parse_args()

    # Get mesh, angles, and parameters
    input_dir = args.input_dir
    mesh_name = args.fname
    output_dir = args.output_dir
    mesh_filepath = os.path.join(input_dir, mesh_name + '.obj')
    v, tc, n, f, ftc, fn = igl.read_obj(mesh_filepath)
    cones = np.loadtxt(os.path.join(input_dir, mesh_name + "_cones"), dtype=int)
    Th_hat = np.loadtxt(os.path.join(input_dir, mesh_name + "_Th_hat"), dtype=float)

    # Create output directory for the mesh
    output_dir = os.path.join(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Load mesh
    ms.clear()
    ms.load_new_mesh(mesh_filepath)

    # Just copy meshes without cones
    if (len(cones) == 0):
        output_path = os.path.join(output_dir, mesh_name + ".obj")
        ms.save_current_mesh(output_path)
        output_angles_path = os.path.join(output_dir, mesh_name+'_Th_hat')
        np.savetxt(output_angles_path, Th_hat)
    else:
        # Select cone vertices
        condselect = "((vi==" + str(cones[0]) + ")"
        for cone in cones[1:]:
            condselect = condselect + " || (vi==" + str(cone) + ")"
        condselect = condselect + ")"

        # Apply meshlab filters
        ms.apply_filter("compute_selection_by_condition_per_vertex", condselect=condselect)
        ms.apply_filter("compute_selection_transfer_vertex_to_face", inclusive=False) 
        ms.apply_filter("apply_selection_dilatation")
        ms.apply_filter("meshing_surface_subdivision_loop", iterations=1, loopweight=0, selected=True)

        # Save mesh to file
        output_path = os.path.join(output_dir, mesh_name + ".obj")
        ms.save_current_mesh(output_path)

        # Get new cone angles for parameterization by adding flat vertex constraints for new vertices
        vr, _, _, _, _, _ = igl.read_obj(output_path)
        num_new_vertices = len(vr) - len(v)
        Th_hat_refined = np.concatenate((Th_hat, np.full(num_new_vertices, 2 * np.pi)))

        # Save final angle to file
        output_angles_path = os.path.join(output_dir, mesh_name+'_Th_hat')
        np.savetxt(output_angles_path, Th_hat_refined)
