# First we install scipy...
import bpy, sys, ensurepip, subprocess, os
ensurepip.bootstrap()  # installs pip into Blender's Python environment
subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "pip"])
user_scripts = bpy.utils.script_path_user()
modules_dir = os.path.join(user_scripts, "modules")
os.makedirs(modules_dir, exist_ok=True)
subprocess.check_call([sys.executable, "-m", "pip", "install",
                       "--target", modules_dir,
                       "scipy"])  # replace with your package
sys.path.append(modules_dir)

from numpy import array, ndindex, ceil, clip
import numpy as np
import importlib
import bpy

# Change script dir below to the path to the PySrc directory
script_dir = "/Users/janba/CloudSrc/DivFree-VectorNoise/PySrc"

# Add the script directory to sys.path if it's not already there
if script_dir not in sys.path:
    sys.path.append(script_dir)

import dfvn
importlib.reload(dfvn)
DFVN_trace = dfvn.DFVN_trace

plo, phi = array([-30, -30, 0]), array([30,  30, 60])

dfvn = DFVN_trace(seed=42, dimensions=3, scale=0.05)
do_reproject = False # Change to true to do reprojection as discussed in paper
N=300
step_len=1.0

# --- Create a Blender object which is a collection of curves from traced paths ---
if 'bpy' in globals():
    # Collect all traced paths and their starting points for color
    paths = []
    for i in range(100):
        p_start = np.array(np.random.uniform(plo, phi))
        path = dfvn.dfvn_multi_trace(p_start, _t=step_len, N=N, w_project=do_reproject)
        paths.append(path)

    # Print bounding box of all paths and check for NaN values
    all_points = np.concatenate(paths, axis=0)
    min_bb = all_points.min(axis=0)
    max_bb = all_points.max(axis=0)
    print("Paths bounding box:", min_bb, max_bb)
    if np.isnan(all_points).any():
        print("Warning: NaN values detected in path points!")

    # Create a new curve data block
    curve_data = bpy.data.curves.new('TracedPaths', type='CURVE')
    curve_data.dimensions = '3D'

    # Create a material for each curve with color based on its starting point
    for i, path in enumerate(paths):
        spline = curve_data.splines.new('POLY')
        spline.points.add(len(path)-1)
        for j, pt in enumerate(path):
            spline.points[j].co = (pt[0], pt[1], pt[2], 1)

    curve_obj = bpy.data.objects.new('streamline_reproject_'+str(do_reproject), curve_data)
    mod = curve_obj.modifiers.new(name="gn", type='NODES')
    mod.node_group = bpy.data.node_groups["make_tubes"]
    bpy.context.collection.objects.link(curve_obj)
