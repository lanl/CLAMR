#!/bin/sh

csplit state.cpp /State::calc_finite_difference/ {5}
mv xx01 calc_finite_difference.cpp
mv xx02 calc_finite_difference_cell_in_place.cpp
mv xx03 calc_finite_difference_face_in_place.cpp
#mv xx04 calc_finite_difference_via_faces.cpp
mv xx05 calc_finite_difference_regular_cells.cpp
#mv xx06 calc_finite_difference_regular_cells_by_faces.cpp
csplit -f yy xx04 /HXRGFLUXIC/
mv yy00 calc_finite_difference_via_faces.cpp
#csplit -f yy xx05 /\*\*\*\*\*\*/
#mv yy00 calc_finite_difference_regular_cells.cpp
csplit -f yy xx06 /HAVE_OPENCL/
mv yy00 calc_finite_difference_regular_cells_by_faces.cpp

cloc --by-file calc*.cpp
lizard calc*.cpp

fgrep -w if calc_finite_difference.cpp |wc
fgrep -w if calc_finite_difference_via_faces.cpp |wc
fgrep -w if calc_finite_difference_call_in_place.cpp |wc
fgrep -w if calc_finite_difference_face_in_place.cpp |wc
fgrep -w if calc_finite_difference_regular_cells.cpp |wc
fgrep -w if calc_finite_difference_regular_cells_by_faces.cpp |wc
