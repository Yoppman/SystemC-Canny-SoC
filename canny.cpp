#include <systemc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265356789202346
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

// Define USE_DEFAULT_FRAMES to use default engineering frames, otherwise use custom frames.
#define USE_DEFAULT_FRAMES  // Comment this line to switch to custom video frames

#ifdef USE_DEFAULT_FRAMES
// Using Default engineering frames
#define COLS 2704
#define ROWS 1520
#define SIZE COLS*ROWS
#define VIDEO_DIR "video_frames"
#define FILENAME "Engineering%03d.pgm"
#else
// Using My video frames
#define COLS 1080
#define ROWS 1920
#define SIZE COLS*ROWS
#define VIDEO_DIR "my_video"
#define FILENAME "Myframe%03d.pgm"
#endif

/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
#define WINSIZE 21

// IMAGE struct for ECPS 203 A5
typedef struct Image_s
{
    unsigned char img[SIZE];

    Image_s(void)
    {
       for (int i=0; i<SIZE; i++)
       { 
          img[i] = 0;
       }
    }

    Image_s& operator=(const Image_s& copy)
    {
       for (int i=0; i<SIZE; i++)
       { 
          img[i] = copy.img[i];
       }            
       return *this;
    }

    operator unsigned char*()
    {
       return img;
    }

    unsigned char& operator[](const int index)
    {
       return img[index];
    }
} IMAGE;

// Utility class for image I/O
class ImageIO
{
public:
    static int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
    {
        FILE *fp;
        char buf[71];
        int r, c;

        if(infilename == NULL) fp = stdin;
        else{
            if((fp = fopen(infilename, "rb")) == NULL){
                fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
                    infilename);
                return(0);
            }
        }

        fgets(buf, 70, fp);
        if(strncmp(buf,"P5",2) != 0){
            fprintf(stderr, "The file %s is not in PGM format in ", infilename);
            fprintf(stderr, "read_pgm_image().\n");
            if(fp != stdin) fclose(fp);
            return(0);
        }
        do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
        sscanf(buf, "%d %d", &c, &r);
        if(c != cols || r != rows){
            fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
            fprintf(stderr, "read_pgm_image().\n");
            if(fp != stdin) fclose(fp);
            return(0);
        }
        do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

        if((unsigned)rows != fread(image, cols, rows, fp)){
            fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
            if(fp != stdin) fclose(fp);
            return(0);
        }

        if(fp != stdin) fclose(fp);
        return(1);
    }

    static int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
        int cols, const char *comment, int maxval)
    {
        FILE *fp;

        if(outfilename == NULL) fp = stdout;
        else{
            if((fp = fopen(outfilename, "wb")) == NULL){
                fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
                    outfilename);
                return(0);
            }
        }

        fprintf(fp, "P5\n%d %d\n", cols, rows);
        if(comment != NULL)
            if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
        fprintf(fp, "%d\n", maxval);

        if((unsigned)rows != fwrite(image, cols, rows, fp)){
            fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
            if(fp != stdout) fclose(fp);
            return(0);
        }

        if(fp != stdout) fclose(fp);
        return(1);
    }
};

// Stimulus module
SC_MODULE(Stimulus)
{
    sc_fifo_out<IMAGE> out_port;

    SC_CTOR(Stimulus)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }

    void process()
    {
        DIR *dir;
        struct dirent *ent;
        char infilename[512];

        // Open video directory
        dir = opendir(VIDEO_DIR);
        if (dir == NULL) {
            fprintf(stderr, "Error opening %s directory.\n", VIDEO_DIR);
            exit(1);
        }

        // Iterate through all files in the video directory
        while ((ent = readdir(dir)) != NULL) {
            // Skip non-pgm files and edge.pgm
            if (strstr(ent->d_name, ".pgm") == NULL || strstr(ent->d_name, "edge.pgm") != NULL) continue;

            // Extract the frame number based on the active filename format
            int file_num;
            if (sscanf(ent->d_name, FILENAME, &file_num) != 1) {
                printf("Filename doesn't match the expected pattern: %s\n", ent->d_name);
                continue; // Skip if the filename doesn't match the expected pattern
            }

            // Use FILENAME and file_num to build the full path
            snprintf(infilename, sizeof(infilename), "%s/" FILENAME, VIDEO_DIR, file_num);

            if(VERBOSE) printf("Reading the image %s.\n", infilename);

            IMAGE image;
            if(ImageIO::read_pgm_image(infilename, image.img, ROWS, COLS) == 0){
                fprintf(stderr, "Error reading the input image, %s.\n", infilename);
                continue;
            }

            // Send the image via out_port
            out_port.write(image);
        }

        // Close the directory
        closedir(dir);
    }
};

// Monitor module
SC_MODULE(Monitor)
{
    sc_fifo_in<IMAGE> in_port;

    SC_CTOR(Monitor)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }

    void process()
    {
        while (true)
        {
            IMAGE edge_image;
            in_port.read(edge_image);

            // Build the output filename
            static int frame_num = 1;
            char outfilename[512];
            snprintf(outfilename, sizeof(outfilename), "%s/%s%03d_edge.pgm", VIDEO_DIR, "Engineering", frame_num++);

            if(VERBOSE) printf("Writing the edge image in the file %s.\n", outfilename);

            if(ImageIO::write_pgm_image(outfilename, edge_image.img, ROWS, COLS, "", 255) == 0){
                fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
                continue;
            }
        }
    }
};

// DataIn module
SC_MODULE(DataIn)
{
    sc_fifo_in<IMAGE> in_port;
    sc_fifo_out<IMAGE> out_port;

    SC_CTOR(DataIn)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }

    void process()
    {
        while (true)
        {
            IMAGE image;
            in_port.read(image);
            out_port.write(image);
        }
    }
};

// DataOut module
SC_MODULE(DataOut)
{
    sc_fifo_in<IMAGE> in_port;
    sc_fifo_out<IMAGE> out_port;

    SC_CTOR(DataOut)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }

    void process()
    {
        while (true)
        {
            IMAGE image;
            in_port.read(image);
            out_port.write(image);
        }
    }
};

// DUT module
SC_MODULE(DUT)
{
    sc_fifo_in<IMAGE> in_port;
    sc_fifo_out<IMAGE> out_port;

    SC_CTOR(DUT)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }

    void process()
    {
        while (true)
        {
            IMAGE image;
            in_port.read(image);

            IMAGE edge;

            // Call canny function
            canny(image.img, ROWS, COLS, SIGMA, TLOW, THIGH, edge.img);

            out_port.write(edge);
        }
    }

    void canny(unsigned char *image, int rows, int cols, float sigma,
             float tlow, float thigh, unsigned char *edge)
    {
        unsigned char nms[SIZE];
        short int smoothedim[SIZE];
        short int delta_x[SIZE];
        short int delta_y[SIZE];
        short int magnitude[SIZE];

        // Initialize arrays
        memset(nms, 0, sizeof(nms));
        memset(smoothedim, 0, sizeof(smoothedim));
        memset(delta_x, 0, sizeof(delta_x));
        memset(delta_y, 0, sizeof(delta_y));
        memset(magnitude, 0, sizeof(magnitude));

        if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
        gaussian_smooth(image, rows, cols, sigma, smoothedim);

        if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
        derivative_x_y(smoothedim, rows, cols, delta_x, delta_y);

        if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
        magnitude_x_y(delta_x, delta_y, rows, cols, magnitude);

        if(VERBOSE) printf("Doing the non-maximal suppression.\n");
        non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);

        if(VERBOSE) printf("Doing hysteresis thresholding.\n");
        apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, edge);
    }

   void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
         short int *magnitude)
   {
      int r, c, pos, sq1, sq2;

      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
            sq1 = (int)delta_x[pos] * (int)delta_x[pos];
            sq2 = (int)delta_y[pos] * (int)delta_y[pos];
            magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
         }
      }
   }

   void derivative_x_y(short int *smoothedim, int rows, int cols,
         short int *delta_x, short int *delta_y)
   {
      int r, c, pos;

      /****************************************************************************
      * Compute the x-derivative. Adjust the derivative at the borders to avoid
      * losing pixels.
      ****************************************************************************/
      if(VERBOSE) printf("   Computing the X-direction derivative.\n");
      for(r=0;r<rows;r++){
         pos = r * cols;
         delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
         pos++;
         for(c=1;c<(cols-1);c++,pos++){
            delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
         }
         delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
      }

      /****************************************************************************
      * Compute the y-derivative. Adjust the derivative at the borders to avoid
      * losing pixels.
      ****************************************************************************/
      if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
      for(c=0;c<cols;c++){
         pos = c;
         delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
         pos += cols;
         for(r=1;r<(rows-1);r++,pos+=cols){
            delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
         }
         delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
      }
   }

   void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
         short int *smoothedim)
   {
      int r, c, rr, cc,     /* Counter variables. */
         windowsize,        /* Dimension of the gaussian kernel. */
         center;            /* Half of the windowsize. */
      float tempim[SIZE]    /* Buffer for separable filter gaussian smoothing. */
         = {0.0},
            kernel[WINSIZE] /* A one dimensional gaussian kernel. */
         = {0.0},
            dot,            /* Dot product summing variable. */
            sum;            /* Sum of the kernel weights variable. */

      /****************************************************************************
      * Create a 1-dimensional gaussian smoothing kernel.
      ****************************************************************************/
      if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
      make_gaussian_kernel(sigma, kernel, &windowsize);
      center = windowsize / 2;

      /****************************************************************************
      * Blur in the x - direction.
      ****************************************************************************/
      if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
      for(r=0;r<rows;r++){
         for(c=0;c<cols;c++){
            dot = 0.0;
            sum = 0.0;
            for(cc=(-center);cc<=center;cc++){
               if(((c+cc) >= 0) && ((c+cc) < cols)){
                  dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
                  sum += kernel[center+cc];
               }
            }
            tempim[r*cols+c] = dot/sum;
         }
      }

      /****************************************************************************
      * Blur in the y - direction.
      ****************************************************************************/
      if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
      for(c=0;c<cols;c++){
         for(r=0;r<rows;r++){
            sum = 0.0;
            dot = 0.0;
            for(rr=(-center);rr<=center;rr++){
               if(((r+rr) >= 0) && ((r+rr) < rows)){
                  dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
                  sum += kernel[center+rr];
               }
            }
            smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
         }
      }
   }

   void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
   {
      int i, center;
      float x, fx, sum=0.0;

      *windowsize = 1 + 2 * ceil(2.5 * sigma);
      center = (*windowsize) / 2;

      if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);

      for(i=0;i<(*windowsize);i++){
         x = (float)(i - center);
         fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
         kernel[i] = fx;
         sum += fx;
      }

      for(i=0;i<(*windowsize);i++) kernel[i] /= sum;

      if(VERBOSE){
         printf("The filter coefficients are:\n");
         for(i=0;i<(*windowsize);i++)
            printf("kernel[%d] = %f\n", i, kernel[i]);
      }
   }

   void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
      int cols)
   {
      short *tempmagptr;
      unsigned char *tempmapptr;
      int i;
      int x[8] = {1,1,0,-1,-1,-1,0,1},
         y[8] = {0,1,1,1,0,-1,-1,-1};

      for(i=0;i<8;i++){
         tempmapptr = edgemapptr - y[i]*cols + x[i];
         tempmagptr = edgemagptr - y[i]*cols + x[i];

         if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
            *tempmapptr = (unsigned char) EDGE;
            follow_edges(tempmapptr,tempmagptr, lowval, cols);
         }
      }
   }

   void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
      float tlow, float thigh, unsigned char *edge)
   {
      int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
      short int maximum_mag=0;

      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
      if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
      else edge[pos] = NOEDGE;
         }
      }

      for(r=0,pos=0;r<rows;r++,pos+=cols){
         edge[pos] = NOEDGE;
         edge[pos+cols-1] = NOEDGE;
      }
      pos = (rows-1) * cols;
      for(c=0;c<cols;c++,pos++){
         edge[c] = NOEDGE;
         edge[pos] = NOEDGE;
      }

      /****************************************************************************
      * Compute the histogram of the magnitude image. Then use the histogram to
      * compute hysteresis thresholds.
      ****************************************************************************/
      for(r=0;r<32768;r++) hist[r] = 0;
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
      if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
         }
      }

      /****************************************************************************
      * Compute the number of pixels that passed the nonmaximal suppression.
      ****************************************************************************/
      for(r=1,numedges=0;r<32768;r++){
         if(hist[r] != 0) maximum_mag = r;
         numedges += hist[r];
      }

      highcount = (int)(numedges * thigh + 0.5);

      r = 1;
      numedges = hist[1];
      while((r<(maximum_mag-1)) && (numedges < highcount)){
         r++;
         numedges += hist[r];
      }
      highthreshold = r;
      lowthreshold = (int)(highthreshold * tlow + 0.5);

      if(VERBOSE){
         printf("The input low and high fractions of %f and %f computed to\n",
      tlow, thigh);
         printf("magnitude of the gradient threshold values of: %d %d\n",
      lowthreshold, highthreshold);
      }

      /****************************************************************************
      * This loop looks for pixels above the highthreshold to locate edges and
      * then calls follow_edges to continue the edge.
      ****************************************************************************/
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
      if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
               edge[pos] = EDGE;
               follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
      }
         }
      }

      /****************************************************************************
      * Set all the remaining possible edges to non-edges.
      ****************************************************************************/
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
      }
   }

   void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
      unsigned char *result)
   {
      int rowcount, colcount,count;
      short *magrowptr,*magptr;
      short *gxrowptr,*gxptr;
      short *gyrowptr,*gyptr,z1,z2;
      short m00; short gx=0; short gy=0;
      float mag1, mag2; float xperp=0;float yperp=0;
      unsigned char *resultrowptr, *resultptr;

      /****************************************************************************
      * Zero the edges of the result image.
      ****************************************************************************/
      for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
         count<ncols; resultptr++,resultrowptr++,count++){
         *resultrowptr = *resultptr = (unsigned char) 0;
      }

      for(count=0,resultptr=result,resultrowptr=result+ncols-1;
         count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
         *resultptr = *resultrowptr = (unsigned char) 0;
      }

      /****************************************************************************
      * Suppress non-maximum points.
      ****************************************************************************/
      for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
         gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
         rowcount<nrows-1; /* bug fix 10/05/23, RD */
         rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
         resultrowptr+=ncols){
         for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
            resultptr=resultrowptr;colcount<ncols-1; /* bug fix 10/05/23, RD */
            colcount++,magptr++,gxptr++,gyptr++,resultptr++){
            m00 = *magptr;
            if(m00 == 0){
               *resultptr = (unsigned char) NOEDGE;
            }
            else{
               xperp = -(gx = *gxptr)/((float)m00);
               yperp = (gy = *gyptr)/((float)m00);
            }

            if(gx >= 0){
               if(gy >= 0){
                     if (gx >= gy)
                     {
                           /* 111 */
                           /* Left point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr - ncols - 1);

                           mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr + ncols + 1);

                           mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                     }
                     else
                     {
                           /* 110 */
                           /* Left point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols - 1);

                           mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                           /* Right point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols + 1);

                           mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                     }
                  }
                  else
                  {
                     if (gx >= -gy)
                     {
                           /* 101 */
                           /* Left point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr + ncols - 1);

                           mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                           /* Right point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr - ncols + 1);

                           mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                     }
                     else
                     {
                           /* 100 */
                           /* Left point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols - 1);

                           mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols + 1);

                           mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                     }
                  }
               }
               else
               {
                  if ((gy = *gyptr) >= 0)
                  {
                     if (-gx >= gy)
                     {
                           /* 011 */
                           /* Left point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr - ncols + 1);

                           mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr + ncols - 1);

                           mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                     }
                     else
                     {
                           /* 010 */
                           /* Left point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols + 1);

                           mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                           /* Right point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols - 1);

                           mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                     }
                  }
                  else
                  {
                     if (-gx > -gy)
                     {
                           /* 001 */
                           /* Left point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr + ncols + 1);

                           mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                           /* Right point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr - ncols - 1);

                           mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                     }
                     else
                     {
                           /* 000 */
                           /* Left point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols + 1);

                           mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols - 1);

                           mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                     }
                  }
               }

               /* Now determine if the current point is a maximum point */

               if ((mag1 > 0.0) || (mag2 > 0.0))
               {
                  *resultptr = (unsigned char) NOEDGE;
               }
               else
               {
                  if (mag2 == 0.0)
                     *resultptr = (unsigned char) NOEDGE;
                  else
                     *resultptr = (unsigned char) POSSIBLE_EDGE;
               }
         }
      }
   }
};

// Platform module
SC_MODULE(Platform) {
    sc_fifo_in<IMAGE> in_port;
    sc_fifo_out<IMAGE> out_port;

    // Change from pointers to instances
    DataIn din;
    DataOut dout;
    DUT duty;
    
    // Change from pointers to instances
    sc_fifo<IMAGE> q1;
    sc_fifo<IMAGE> q2;

    SC_CTOR(Platform)
    : din("din")
    , duty("duty")
    , dout("dout")
    , q1("q1_internal", 1)
    , q2("q2_internal", 1)
    {
        // Bind ports using references
        din.in_port.bind(in_port);
        din.out_port.bind(q1);

        duty.in_port.bind(q1);
        duty.out_port.bind(q2);

        dout.in_port.bind(q2);
        dout.out_port.bind(out_port);
    }
};

// Top module
SC_MODULE(Top) {
    // Change from pointers to instances
    Stimulus stimulus;
    Platform platform;
    Monitor monitor;
    
    // Change from pointers to instances
    sc_fifo<IMAGE> q1;
    sc_fifo<IMAGE> q2;

    SC_CTOR(Top)
    : stimulus("stimulus")
    , platform("platform")
    , monitor("monitor")
    , q1("q1", 1)
    , q2("q2", 1)
    {
        // Bind ports using references
        stimulus.out_port.bind(q1);
        platform.in_port.bind(q1);
        platform.out_port.bind(q2);
        monitor.in_port.bind(q2);
    }
};

int sc_main(int argc, char* argv[])
{
    Top top("top");

    sc_start();

    return 0;
}