#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>
#include <stdint.h>
#include <float.h>


/*****************************************************************************
 * 	kahan_sum()
 * 	Add a list of values using the Kahan method
 *
 * 	Params:
 * 	f_list - list of real_ts to sort
 * 	size - size of the list
 *****************************************************************************/
inline real_t kahan_sum(real_t *f_list, int size){
    real_t sum = 0.0;
    real_t c = 0.0;
    real_t t, y;
    int i;
    for (i = 0; i < size; i++){
        y = f_list[i] - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}


/*****************************************************************************
 * 	psum_add()
 *      Add a list of values using the Kahan method
 *
 * 	Params:
 * 	f_list - list of real_ts to sort
 * 	size - size of the list
 ******************************************************************************/
inline real_t psum(real_t *f_list, int size){
    real_t sum = 0.0;
#ifdef FULL_PRECISION
    real_t smallest = DBL_MAX;
#else
    real_t smallest = FLT_MAX;
#endif
    real_t temp_sum1, temp_sum2;

    int left = -1;
    int right = -1;
    int center = -1;
    int i;

    for (i = 0; i < size; i++){
        if (fabs(f_list[i]) < fabs(smallest)){
            smallest = f_list[i];
            center = i;
        }
    }

    left = center -1;
    right = center + 1;
    sum = f_list[center];

    for (i = 0; i < size - 1; i++){

        // Hit the end of left, add down right
        if (left < 0){
           sum += f_list[right];
           right++;
        } 
        // Hit the end of right, add down left
        else if (right > size){
           sum += f_list[left];
           left--;
        }
        // Compare both, add smallest
        else {
            temp_sum1 = sum + f_list[left];
            temp_sum2 = sum + f_list[right];

            fabs(temp_sum1) < fabs(temp_sum2) ? left-- : right++;
            sum = (fabs(temp_sum1) < fabs(temp_sum2)) ? temp_sum1 : temp_sum2;
        }
    }
    return sum;
}


/*****************************************************************************
	sort_list()
	Sorts a list of real_ts based on order of magnitude
	
	Params:
		f_list - list of real_ts to sort
		size - size of the list
		order - 'a': ascending will sort from small magnitude to larger magnitude
		      - 'd': descending will sort from large magnitude to smaller magnitude
*****************************************************************************/
inline void sort_list(real_t *f_list, int size, char order){
	int i, j;
	real_t temp;
	int sign_change = 1;

	if (order == 'd')
		sign_change = -1;

	for(i = 0; i < size; i++){
		for(j = 0; j < size-i-1; j++){
			if (f_list[j] * sign_change > f_list[j+1] * sign_change ){
				temp = f_list[j+1];
				f_list[j+1] = f_list[j];
				f_list[j] = temp;
			}
		}
	}
}

inline real_t U_fullstep_version(
        real_t    deltaT,
        real_t    dr,
        real_t    U,
        real_t    F_plus,
        real_t    F_minus,
        real_t    G_plus,
        real_t    G_minus,
        real_t    wplusx_H,
        real_t    wminusx_H,
        real_t    wplusy_H,
        real_t    wminusy_H
) {

real_t f_list[9] = {U, -(deltaT/dr)*F_plus, -(deltaT/dr)*(-F_minus), -(deltaT/dr)*(G_plus),-(deltaT/dr)*(-G_minus), -wminusx_H, wplusx_H, -wminusy_H, wplusy_H};

#ifdef UPDATE_EQUATION_SORT_DESCENDING
sort_list(f_list,9,'d');
#endif

#ifdef UPDATE_EQUATION_SORT_ASCENDING
sort_list(f_list,9,'a');
#endif

#ifdef UPDATE_EQUATION_KAHAN
return kahan_sum(f_list,9);
#endif

#ifdef UPDATE_EQUATION_PSUM
return psum(f_list,9);
#endif

#ifdef UPDATE_EQUATION_PAIRWISE
return (((U+wplusx_H)+((deltaT/dr)*F_minus+-wminusx_H))+(((-deltaT/dr)*F_plus+-wminusy_H)+((-deltaT/dr)*G_plus+((deltaT/dr)*G_minus+wplusy_H))));
#endif
}
