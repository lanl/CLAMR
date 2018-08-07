/*
 *  Copyright (c) 2011-2012, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Copyright 2011-2012. Los Alamos National Security, LLC. This software was produced 
 *  under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
 *  ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
 *  ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 *  to produce derivative works, such modified software should be clearly marked,
 *  so as not to confuse it with the version available from LANL.
 *
 *  Additionally, redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE LOS ALAMOS NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
 *  SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *  
 *  CLAMR -- LA-CC-11-094
 *  This research code is being developed as part of the 
 *  2011 X Division Summer Workshop for the express purpose
 *  of a collaborative code for development of ideas in
 *  the implementation of AMR codes for Exascale platforms
 *  
 *  AMR implementation of the Wave code previously developed
 *  as a demonstration code for regular grids on Exascale platforms
 *  as part of the Supercomputing Challenge and Los Alamos 
 *  National Laboratory
 *  
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Daniel Dunning  XCP-2   ddunning@lanl.gov, daniel.dunning@ttu.edu
 * 
 */

#define is_lower(i) (i % 2 == 0) 
#define is_upper(i) (i % 2 == 1) 

__kernel void calc_face_list_wbidirmap(
                        const int       ncells,                 // 0
                        const int       nxfaces,                // 1
                        const int       nyfaces,                // 2
                        const int       levmx,                  // 3
            __global          int   *level,                     // 4
            __global          int   *nlft,                      // 5
            __global          int   *nrht,                      // 6
            __global          int   *nbot,                      // 7
            __global          int   *ntop,                      // 8
            __global          int   *map_xface2cell_lower,      // 9
            __global          int   *map_xface2cell_upper,      // 10
            __global          int   *map_xcell2face_left1,      // 11
            __global          int   *map_xcell2face_left2,      // 12
            __global          int   *map_xcell2face_right1,     // 13
            __global          int   *map_xcell2face_right2,     // 14
            __global          int   *xface_level,               // 15
            __global          int   *xface_i,                   // 16
            __global          int   *xface_j,                   // 17
            __global          int   *ixmin_level,               // 18
            __global          int   *ixmax_level,               // 19
            __global          int   *jxmin_level,               // 20
            __global          int   *jxmax_level,               // 21
            __global          int   *map_yface2cell_lower,      // 22
            __global          int   *map_yface2cell_upper,      // 23
            __global          int   *map_ycell2face_bot1,       // 24
            __global          int   *map_ycell2face_bot2,       // 25
            __global          int   *map_ycell2face_top1,       // 26
            __global          int   *map_ycell2face_top2,       // 27
            __global          int   *yface_level,               // 28
            __global          int   *yface_i,                   // 29
            __global          int   *yface_j,                   // 30
            __global          int   *iymin_level,               // 31
            __global          int   *iymax_level,               // 32
            __global          int   *jymin_level,               // 33
            __global          int   *jymax_level) {             // 34

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const unit giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const unit ntX = get_local_sizee(0);

    const unit group_id = get_group_id(0);

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    if (giX >= max(ncells, max(nxface, nyface))
        return;

    // Ensure the executing thread is not extraneous
    if (giX < ncells) {
        int iface = xfaceIdxList[giX], nz = giX;
        int nr = nrht[nz];
        if (nr != nz) {

            int ifactor = 1;
            if (level[nr] < level[nz]) ifactor = 2;

            map_xface2cell_lower[iface] = nz;
            map_xface2cell_upper[iface] = nr;
            xface_levl[iface] = MAX(level[nz], level[nr]);
            xface_i[iface] = i[nr]*ifactor;
            if (level[nr] < level[nz] && is_upper(j[nz]) ) {
                xface_j[iface] = j[nr]*ifactor+1;
            }
            else {
                xface_j[iface] = j[nr]*ifactor;
            }
            map_xcell2face_right1[nz] = iface;

            iface++; 

            if (level[nr] > level[nz] && is_lower(j[nr]) ){
                int ntr = ntop[nr];
                if (ntr != nr) {
                    iface++;
                    map_xface2cell_lower[iface] = nz;
                    map_xface2cell_upper[iface] = ntr;
                    xface_level[iface] = MAX(level[nz],level[ntr]);
                    xface_i[iface] = i[ntr]*ifactor;
                    xface_j[iface] = j[ntr]*ifactor;
                    map_xcell2face_right2[nz] = iface;

                    iface++;
                }
            }
        }

        int nl = nlft[nz];
        if (nl != nz) { 
            if (level[nl] < level[nz] && is_upper(j[nz])){
                 map_xcell2face_left1[nz] = map_xcell2face_right2[nl];
            } else {
                map_xcell2face_left1[nz] = map_xcell2face_right1[nl];
                if (level[nl] > level[nz]){
                    map_xcell2face_left2[nz] = map_xcell2face_right1[ntop[nl]];
                }
            }
        }

        //for yfaces
        iface = yfaceIdxList[giX], nz = giX;
        int nt = ntop[nz];
        if (nt != nz) {

            int ifactor = 1;
            if (level[nt] < level[nz]) ifactor = 2;

            map_yface2cell_lower[iface] = nz;
            map_yface2cell_upper[iface] = nt;
            yface_levl[iface] = MAX(level[nz], level[nt]);
            yface_i[iface] = j[nt]*ifactor;
            if (level[nt] > level[nz] && is_upper(i[nz]) ) {
                yface_i[iface] = i[nt]*ifactor+1;
            }
            else {
                yface_i[iface] = i[nt]*ifactor;
            }
            map_ycell2face_top1[nz] = iface;

            iface++; 

            if (level[nt] > level[nz] && is_lower(i[nt]) ){
                int nrt = nrht[nt];
                if (nrt != nt) {
                    iface++;
                    map_yface2cell_lower[iface] = nz;
                    map_yface2cell_upper[iface] = nrt;
                    yface_level[iface] = MAX(level[nz],level[nrt]);
                    yface_i[iface] = i[nrt]*ifactor;
                    yface_j[iface] = j[nrt]*ifactor;
                    map_ycell2face_top2[nz] = iface;

                    iface++;
                }
            }
        }

        int nb = nbot[nz];
        if (nb != nz) { 
            if (level[nb] < level[nz] && is_upper(i[nz])){
                 map_ycell2face_bot1[nz] = map_ycell2face_top2[nb];
            } else {
                map_xcell2face_bot1[nz] = map_xcell2face_top1[nb];
                if (level[nb] > level[nz]){
                    map_ycell2face_bot2[nz] = map_ycell2face_top1[nrht[nb]];
                }
            }
        }
    }
        
    if (giX < nxface) {
        int fl = xface_level[iface];

        int fi = xface_i[iface];
        if (fi < ixmin_level[fl]) ixmin_level[fl] = fi;
        if (fi > ixmax_level[fl]) ixmax_level[fl] = fi;

        int fj = xface_j[iface];
        if (fj < jxmin_level[fl]) jxmin_level[fl] = fj;
        if (fj > jxmax_level[fl]) jxmax_level[fl] = fj;
    
        int fl = xface_level[iface];
        if (ixmax_level[fl] > ixmin_level[fl]) {
            xface_i[iface] -= ixmin_level[fl];
            xface_j[iface] -= jxmin_level[fl];
        }
    }

    if (giX < nyface) {
        int fl = yface_level[iface];

        int fi = yface_i[iface];
        if (fi < iymin_level[fl]) iymin_level[fl] = fi;
        if (fi > iymax_level[fl]) iymax_level[fl] = fi;

        int fj = yface_j[iface];
        if (fj < jymin_level[fl]) jymin_level[fl] = fj;
        if (fj > jymax_level[fl]) jymax_level[fl] = fj;
        
        int fl = yface_level[iface];
        if (iymax_level[fl] > iymin_level[fl]) {
            yface_i[iface] -= iymin_level[fl];
            yface_j[iface] -= jymin_level[fl];
        }
    }

    if (giX < levmx) {
        ixadjust[fl] = ixmin_level[fl];
        jxadjust[fl] = jxmin_level[fl];
        ixmax_level[fl] -= ixmin_level[fl];;
        jxmax_level[fl] -= jxmin_level[fl];
        ixmin_level[fl] = 0;
        jxmin_level[fl] = 0;

        iyadjust[fl] = iymin_level[fl];
        jyadjust[fl] = jymin_level[fl];
        iymax_level[fl] -= iymin_level[fl];;
        jymax_level[fl] -= jymin_level[fl];
        iymin_level[fl] = 0;
        jymin_level[fl] = 0;
    }
}

//in the regular mesh class, will need dev memory for x/yfaceIdxList
/*for (int ccc = 0; ccc < (int) ncells - 1; ccc++) {
    int xfaceSt = xfaceIdxList[ccc], yfaceSt = yfaceIdxList[ccc];
    xfaceSt++;
    nxface++;
    yfaceSt++;
    nyface++;
    if (level[nrht[ccc]] > level[ccc]) {
        xfaceSt++; 
        nxface++;
    }
    if (level[ntop[ccc]] > level[ccc]) {
        yfaceSt++;
        nyface++;
    }
    xfaceIdxList[ccc+1] = xfaceSt;
    yfaceIdxList[ccc+1] = yfaceSt;
}*/


