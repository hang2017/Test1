#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include "post_process.h"

void moy_lgn(double *tab, int width, int height, double *output) {
	int i,j;

	for(i = 0;i < height;i++){
		double res = 0;
		for(j = 0; j<width;j++){
			res += tab[i*width + j];
		}
		res = res/width;
		output[i] = res;
	}
}

void avg_row_parallel(double *tab, int width, int height, double *output, int margin, MPI_Comm comm) {
    /* Computes the row-by-row average in parallel.
     *
     * Assumes that MPI_Comm is a communicator of a cartesian grid of process
     * Computes the local average row-by-row in each process
     * Forwards the result along with the number of elements averaged to the leftmost processes
     * The leftmost processes then computes the global row-by-row average
     */
    int i, j;
    double res;
    int rank;
    MPI_Comm row_comm;
    int *widths;
    double *gathered_values;
    int dims[2], periods[2], car_coord[2];
    int elt_sum = 0;

    // True width and height
    int theight = height-2*margin;
    int twidth = width-2*margin;

    // Get the process grid parameters
    MPI_Cart_get(comm, 2, dims, periods, car_coord);

    assert(margin >= 0);
    // Computes the local average, store the result in output
    for (i = margin; i < height - margin; i++) {
        res = 0;
        for (j = margin; j < width - margin; j++) {
            res += tab[i*width + j];
        }
        res = res/twidth;
        output[i-margin] = res;
    }

    // Split the process grid by line
    MPI_Comm_split(comm, car_coord[1], car_coord[0], &row_comm);
    MPI_Comm_rank(row_comm, &rank);

    // The leftmost process in each line have rank 0
    if (rank == 0) {
        widths = malloc(dims[0]*sizeof(int));
        gathered_values = malloc(theight*dims[0]*sizeof(double));
    }

    // Gather the results of the processes in the same line
    MPI_Gather(output, theight, MPI_DOUBLE, gathered_values, theight, MPI_DOUBLE, 0, row_comm);
    MPI_Gather(&twidth, 1, MPI_INT, widths, 1, MPI_INT, 0, row_comm);

    // Computes global average
    if (rank == 0) {
        for (i=0; i < theight; i++) output[i] = 0; //Clean output since it's dirty for the leftmost processes
        for (j=0; j < dims[0]; j++) {
            elt_sum += widths[j];
            for (i=0; i < theight; i++) {
                output[i] += gathered_values[j*theight+i]*widths[j];
            }
        }
        for (i=0; i < theight; i++) {
            assert(elt_sum > 0);
            output[i] = output[i]/elt_sum;
        }
        free(widths);
        free(gathered_values);
    }
}

void moy_col(double *tab, int width, int height,double *output){
	int i,j;

	for(i = 0;i < width;i++){
		double res=0;
		for(j = 0; j<height;j++){
			double x = tab[j*width + i];
			res += x;
		}
		res  = res/height;
		output[i] = res;
	}

}

double moy_global(double *tab, int width, int height){
	int i,j;
	double res = 0;

	for(i =0; i<height; i++){			//global
		for(j = 0; j<width;j++){
			res += tab[i*width + j];
		}
	}
	res = res/height/width;

	return res;
}

void derive(double *tab, double *next, int width, int height, double *output){
	int i,j;
	double res = 0;

	for(i =0 ;i<height; i++){
		for(j = 0;j<width;j++){
			double t1,t2;
			t1 = next[i*width+j];
			t2 = tab[i*width+j];
			res = t1 - t2;

			output[i*width + j] = res;
		}
	}
}

void laplacien(double *tab, int width, int height,double *output){
	double res = 0;

	for(int i = 1;i < height-1;i++){
		for(int j = 1; j < width-1; j++){
			res = 4*tab[i*width+j]-tab[i*width + j-1]-tab[i*width + j+1]-tab[(i+1)*width + j]-tab[(i-1)*width + j];
			output [i*width + j] = res;
		}
	}
}
