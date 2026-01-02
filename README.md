# DivFree-VectorNoise
This repository contains shaders and Python code that demonstrates methods from the paper "Improving Curl Noise" which is available here: https://doi.org/10.1145/3757377.3763980

## Shaders

The `.fs` files in the `Shaders` directory are meant to be run in ShaderToy.

- `image_warping.fs` produces the warped images in Figure 5.

## Python and Blender files

In the `PySrc` directory you will find a collection of Python files and a blender file (for Blender 5.0.0). If you open the Blender file (`streamlines.blend`) you should see two collections of stream lines (streamline_reproject_False and streamline_reproject_True) in the object hierarchy. The object names indicate whether our **reprojection** scheme was applied. When reproject is false the streamlines do not close up and when it is true many of the streamlines do close up. The streamlines are always closed, but since we only trace a finite length, the segment that is traced may not come full circle.

You can compute the streamlines yourself by loading and running the script `streamlines.py` from within Blender. The script requires `scipy` and should install it on its own. In line 19 you need to set the path to the `PySrc` directory before running otherwise key python files are not available. In line 32 you indicate whether to do reprojection. Other parameters are also available for tuning. Feel free to experiment.

# Authors

_Andreas Bærentzen, Jonas Martinez, Jeppe R. Frisvad, Sylvain Lefebvre_
