#!/bin/bash

outlines=0
FPS=6
SIZE=800

while getopts ys:f: opt 
   do
      case "$opt" in
         y) outlines=1;;
         s) SIZE="$OPTARG";;
         f) FPS="$OPTARG";;
         [?]) exit 1;;
      esac
done

if [ $outlines -eq 1 ]; then
   python generate_image.py y
else
   python generate_image.py
fi

cd ./graphics_output

echo $SIZE

mencoder mf://@imagelist.txt -mf w=800:h=800:fps=$FPS:type=png -vf scale=$SIZE:$SIZE -oac copy -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:vbitrate=8000000 -o test.avi
#mencoder mf://@imagelist.txt -mf w=$SIZE:h=$SIZE:fps=$FPS:type=png -ovc copy -oac copy -o test.avi
mv test.avi ../
