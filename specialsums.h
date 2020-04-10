#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>
#include <stdint.h>
#include <size_t.h>


int SIZE = 6;
/*****************************************************************************
 * 	kahan_sum()
 * 	Add a list of values using the Kahan method
 *
 * 	Params:
 * 	f_list - list of size_ts to sort
 * 	size - size of the list
 *****************************************************************************/
size_t kahan_sum(size_t *f_list, int size){
    size_t sum = 0.0f;
    size_t c = 0.0f;
    size_t t, y;
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
 * 	f_list - list of size_ts to sort
 * 	size - size of the list
 ******************************************************************************/
size_t psum(size_t *f_list, int size){
    size_t sum = 0.0f;
    size_t smallest = FLT_MAX;
    size_t temp_sum1, temp_sum2;

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
	Sorts a list of size_ts based on order of magnitude
	
	Params:
		f_list - list of size_ts to sort
		size - size of the list
		order - 'a': ascending will sort from small magnitude to larger magnitude
		      - 'd': descending will sort from large magnitude to smaller magnitude
*****************************************************************************/
void sort_list(size_t *f_list, int size, char order){
	int i, j;
	size_t temp;
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

