#!/usr/bin/python
import sys, glob, os
from PIL import Image, ImageDraw, ImageColor, ImageFont

def check_outline_flag(arg):
   if arg == 'y':
      return True
   else:
      print(arg + " is an unrecognized command")
      print("If you want grid lines pass y to generate_image.py")
      sys.exit()

def write_cp_img_to_list(cp_num, img_list_file, num_cp_imgs, final_graph_num):
   total_writes = num_cp_imgs * cp_num
   start_graph = final_graph_num - total_writes
   last_graph = start_graph + total_writes
   if start_graph == 0:
      start_graph += 1
   if cp_num == 1:
      img_list_file.write('image'+str(final_graph_num)+'cp'+str(cp_num)+'.png\n')
      for i in range(start_graph,last_graph):
         img_list_file.write('image'+str(i)+'cp'+str(cp_num)+'.png\n')
      return
   else:
      write_cp_img_to_list(cp_num - 1, img_list_file, num_cp_imgs, final_graph_num)
      img_list_file.write('image'+str(final_graph_num)+'cp'+str(cp_num)+'.png\n')
      for i in range(start_graph,last_graph):
         img_list_file.write('image'+str(i)+'cp'+str(cp_num)+'.png\n')
      return


GridLines = False
add_one = False
NumImages=0
Iteration = 0
Time = 0.0
font = ImageFont.truetype('FreeSansBold.ttf',20)
cp_info = [-1,-1,-1] # [Checkpoint Num, Last Graph for Checkpoint Num, Total number of Checkpoint Images for Checkpoint]

if len(sys.argv) == 2:
   GridLines = check_outline_flag(sys.argv[1])
else:
   GridLines = False
   
os.chdir('./graphics_output')
for file in glob.glob('*.data'):
   im = Image.new("RGB",(800,800))
   draw = ImageDraw.Draw(im)
   name = file
   file_num = name.strip('graph')
   file_num = file_num.strip('.data')
   found = file_num.find('cp')
   if(found > 0):
      graph_num = int(file_num[0:found])
      cp_num = int(file_num[found+2:])
      if cp_info[0] < cp_num :
         cp_info[0] = cp_num
         if cp_info[1] < graph_num:
            cp_info[1] = graph_num
         cp_info[2] = 1
      elif cp_info[0] == cp_num:
         if cp_info[1] < graph_num:
            cp_info[1] = graph_num
         cp_info[2] += 1
      if graph_num == 1:
         add_one = True
   with open(name) as f:
      for line in f:
         line = line.strip('\n')
         cellInfo = line.split(',')
         if len(cellInfo) == 2:
            Iteration = cellInfo[0]
            Time = cellInfo[1]
         else:
            for i in range (len(cellInfo)):
               cellInfo[i] = int(cellInfo[i])
            x1 = cellInfo[0]
            x2 = cellInfo[0] + cellInfo[2]
            y1 = cellInfo[1]
            y2 = cellInfo[1] + cellInfo[3]
            color = str(cellInfo[4])
            fill_str = 'hsl(' + color + ',100%,50%)'
            draw.rectangle([x1,y1,x2,y2],fill=fill_str)
   if GridLines:
      outline_file = 'outline' + file_num + '.lin'
      with open(outline_file) as f:
         for line in f:
            line = line.strip('\n')
            cellInfo = line.split(',')
            for i in range (len(cellInfo)):
               cellInfo[i] = int(cellInfo[i])
            x1 = cellInfo[0]
            x2 = cellInfo[2]
            y1 = cellInfo[1]
            y2 = cellInfo[3]
            draw.line([x1,y1,x2,y2],fill=(0,0,0))
   draw.text([0,20],"Iteration: " + Iteration,font=font,fill=(255,255,255))
   draw.text([0,40],"Time: " + Time,font=font,fill=(255,255,255))
   name = 'image' + str(file_num) +'.png'
   im.save(name,'PNG')
   NumImages+=1
   print name + ' has been generated'

if add_one:
   num_cp_images = cp_info[2]/cp_info[0]
else:
   num_cp_images = (cp_info[2] - 1)/cp_info[0]
print 'The number of cp images are ' + str(num_cp_images) 
with open('imagelist.txt','w') as fo:
   for i in range(0,NumImages):
      if i == cp_info[1]:
         write_cp_img_to_list(cp_info[0],fo,num_cp_images,cp_info[1])
      fo.write('image%05d.png\n' % (i))
