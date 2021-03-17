// -----------------------------------------------------------------
// Overlapping 2D-DCT (Discrete Cosine Transform) image filter
//
// Image model: z = y + n, where z is measured image mxn
//                               y is the noise free image
//                               n is Gaussian noise
//                               
// 
// Compile with: g++ -O2 filter.c dct.c -o filter
//
// Mikael Mieskolainen, 2011
// ----------------------------------------------------------------- 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dct.h" // 2D-DCT Library


// This abstracts the image matrix as linear array
// access image matrix pixel at row i and column j, with N columns
// The first for "->" and the second for "." access
#define pixel(data,i,j)   ( (data)->image[(i)*((data)->N) + (j)] )
#define pixel_h(data,i,j) ( (data).image[(i)*((data).N) + (j)] )


// This is the image object
struct IMAGE {

  IMAGE(int m, int n) {
    M = m;
    N = n;
  };

  double *image; // Data
  int M;         // Number of Rows
  int N;         // Number of Columns
};


// Subfunction prototypes
void ones(struct IMAGE*);
void zeros(struct IMAGE*);
void printValues(struct IMAGE*);
void pointDiv(struct IMAGE*, struct IMAGE*, struct IMAGE*);
void pointMul(struct IMAGE*, struct IMAGE*, struct IMAGE*);
void dctFilter(struct IMAGE*, struct IMAGE*, int hop, double thresh, int window_size);


// Main program
int main(int argc, char* argv[]) {

  // Setup the parameters
  char inputfile[]   = "lena_noisy.bin";
  char outputfile[]  = "lena_denoised.bin";

  // Image dimensions fixed here (no header in the binary files)
  int M = 512; // Rows
  int N = 512; // Cols
  
  
  // Filter parameters
  int hop            = 1;     // Local window hopping step
  double threshold   = 40;    // Frequency domain threshold

  const int window_size = 8;  // Keep it at 8 (2D-DCT function is fixed 8x8)


  // ------------------------------------------------

  // Create image matrix for the original (noisy) image
  // and denoised image
  struct IMAGE z(M,N);
  struct IMAGE y_est(M,N);

  // ------------------------------------------------

  // Allocate memory
  z.image = (double*) malloc(z.M * z.N * sizeof(double));
  y_est.image = (double*) malloc(y_est.M * y_est.N * sizeof(double));
  
  if (z.image == NULL) {
    printf("Out of memory with z \n");
    return EXIT_FAILURE;
  }
  if (y_est.image == NULL) {
    printf("Out of memory with y_est \n");
    return EXIT_FAILURE;
  }

  // Init matrices with zeros
  zeros(&z);
  zeros(&y_est);
  
  // ------------------------------------------------
  // Read the image
  
  FILE *fid;

  fid = fopen(inputfile, "rb");
  if (fid) {
    int ret = fread(z.image, 8, z.M * z.N, fid);
    fclose(fid);
  } else {
    printf("Error opening input image! \n");
    return EXIT_FAILURE;
  }
  
  // ------------------------------------------------
  // Filter the image

  dctFilter(&z, &y_est, hop, threshold, window_size);
  
  // ------------------------------------------------
  // Write the image
  fid = fopen(outputfile, "wb");
  if (fid) { 
    int ret = fwrite(y_est.image, 8, y_est.M * y_est.N, fid);
    fclose(fid);
  } else {
    printf("Error opening file for writing! \n");
    return EXIT_FAILURE;
  }
  
  // -----------------------------------------------
  
  printf("2D-DCT filtering done!\n");
  
  // Free the image matrices
  free(z.image);
  free(y_est.image);
    
  return EXIT_SUCCESS;
}


void dctFilter(struct IMAGE* z, struct IMAGE* y_est, int hop, double thresh, int window_size) {
  
  // Boundary around image matrix
  const int b = window_size / 2;

  // 2D-DCT 8x8 buffers (matrix row by row as linear array)
  double inblock[(window_size*window_size)]  = {0};
  double outblock[(window_size*window_size)] = {0};

  // Allocate the aggregation matrix and weight matrix
  struct IMAGE A(z->M, z->N);
  struct IMAGE W(z->M, z->N);

  A.image = (double*) malloc(A.M * A.N * sizeof(double));
  W.image = (double*) malloc(W.M * W.N * sizeof(double));
  
  if (A.image == NULL) {
    printf("Out of memory with A \n");
    return;
  }
  if (W.image == NULL) {
    printf("Out of memory with W \n");
    return;
  }
  
  // Init matrices with zeros
  zeros(&A);
  zeros(&W);
  
  // Over all rows
  for (int i = b; i < z->M - b + 2; i += hop) {

    // Over all columns
    for (int j = b; j < z->N - b + 2; j += hop) {

      double* VP = &inblock[0]; // Current pointer

      // Accumulate the local window
      for (int m = i-b, k = 0; m < i+b; ++m) {
        for (int n = j-b; n < j+b; ++n) {

          // Check boundary
          *(VP + k) = (m < z->M && n < z->N) ? pixel(z,m,n) : 0.0;
          ++k;
        }
      }

      // 2D-DCT for the block
      fdct(VP, outblock);
      
      // ------------------------------------------------------------
      // Filtering the frequency domain & weight calculation
      //
      // N.B. different strategies for the weights are possible

      double sum  = 0.0;
      double sum2 = 0.0;
      for (int p = 0; p < (window_size*window_size); ++p) {

        const double c = fabs(outblock[p]);

        if (c < thresh) {
          outblock[p] = 0.0;
        } else {
          sum  += c;   // Accumulate coefficient values
          sum2 += c*c;
        }
      }
      const double var = (sum2 - sum*sum/window_size) / (window_size*window_size);
      
      // Weight for this block
      const double W_i = var > 0 ? var : 1.0;
      // ------------------------------------------------------------

      // Inverse 2D-DCT for the block
      idct(outblock, VP);

      // Aggregate the block and weight
      for (int m = i-b, k = 0; m < i+b; ++m) {
        for (int n = j-b; n < j+b; ++n) {

          if (m < z->M && n < z->N) { // Check boundary
            pixel_h(A,m,n) += W_i * (*(VP + k));
            pixel_h(W,m,n) += W_i;
          }
          ++k; // here, outside if {}
        }
      }
    }
  }
  
  // Finally divide the aggregation buffer with the weight buffer
  // to normalize the weighted mean
  for (int i = 0; i < y_est->M; ++i) {
    for (int j = 0; j < y_est->N; ++j) {
      pixel(y_est,i,j) = pixel_h(A,i,j) / pixel_h(W,i,j);
    }  
  }
  
  // Free the image matrices
  free(A.image);
  free(W.image);
}


// Point wise division: y = a ./ b
void pointDiv(struct IMAGE *y, struct IMAGE *a, struct IMAGE *b) {
  
  for (int i = 0; i < y->M; ++i) {
    for (int j = 0; j < y->N; ++j) {
      pixel(y,i,j) = pixel(a,i,j) / pixel(b,i,j);
    }  
  }   
}


// Point wise multiplication: y = a .* b
void pointMul(struct IMAGE *y, struct IMAGE *a, struct IMAGE *b) {
  
  for (int i = 0; i < y->M; ++i) {
    for (int j = 0; j < y->N; ++j) {
      pixel(y,i,j) = pixel(a,i,j) * pixel(b,i,j);
    }
  }
}


void ones(struct IMAGE *x) {
  
  for (int i = 0; i < x->M; ++i) {
    for (int j = 0; j < x->N; ++j) {
      pixel(x,i,j) = 1;
    }
  }
}


void zeros(struct IMAGE *x) {
  
  int count = 1;
  for (int i = 0; i < x->M; ++i) {
    for (int j = 0; j < x->N; ++j) {
      pixel(x,i,j) = 0;
    }
  }
}


void printValues(struct IMAGE *x) {
  
  for (int i = 0; i < x->M; i = i + 10) {
    for (int j = 0; j < x->N; j = j + 10) {
      printf("Value at [%d,%d] is %0.2f \n", i, j, pixel(x,i,j) );
    }
  }
}

