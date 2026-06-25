# Spectre-Tile-Creation-Tools
Two Python scripts to create graphic images and printable STL files of Spectre tiles with wiggly edges

This project has two Python scripts:

* spectre_curved_edge_generator.py : This creates graphic images of Spectre tiles, optionally with curved/wiggly edges (is there a convention for describing this?) that when 3D printed allow them to be interlocked, and that create certain merging restrictions in large assemblies of tiles. See the images subdirectory for examples of this script's output.
* spectre_STL_tile_generator.py : This script accepts the user's dimensions and wiggle-depth specifications and creates printable STL files of user-specified sizes. This script also shrinks its tiles to allow them to interlock without too much friction. The shrink amount is under the user's control.

See the example graphics in the images directory, and the STL files with accompanying images in the spectre_tiles directory.

The Spectre tile is a relatively recent creation. It can tile the plane, but without repetition (this is called "aperiodic"). This follows many efforts to create such a tile, beginning with schemes requiring multiple tiles, then a single tile that could tile the plane but had to be flipped over in certain places (a "reflection"), now the Spectre, which does away with that limitation. The Spectre tile is quite an achievement, and just im time to be 3D printed.