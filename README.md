#Material Point Method for Snow
This project implements the material point method for snow as described by Stomakhin et al. 2013.
The code will run a snowball smash between a static larger and a smaller snowball with initial velocity.
Obj files can be loaded in to be the show. There is also a function to create an evenly distributed pointlcoud for a sphere in the code.
The pointcloud (.ptc) files that are created by the program can be converted to bgeo files with the python script attached.
Those can then be loaded into the Houdini example scene and rendered out.
The NGL library, EIGEN library, Renderman pointcloud, and openmp are needed to compile this project.
