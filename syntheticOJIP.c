/*=========================================================================
 * syntheticOJIP.c - generating and binning synthetic photons
 *
 * This routine resamples histograms of photon arrival times
 * It takes the measured bin widths of actual TDC histograms.
 * It simulates randomly distributed photons within these histograms.
 * It then creates a new histogram with equally sized bins.
 *
 * Inputs: inputHist, firstCalBin, lastCalBin, realBinWidth, linearBinWidth
 * Output: corrHist
 *
 * dataStream     - Raw TDC data stream array (MxN uint16 array),
 *                  where M is number of bins, N is number of frames
 * firstCalBin    - First indices of calibrated bins (1xN int16 vector)
 *                  Calibrated bins are those in the TDC range containing 
 *                  useful data
 * lastCalBin     - Last indices of calibrated bins (1xN int16 vector)
 * realBinWidth   - Array of calibrated bin widths (MxN double array)
 * linearBinWidth - Linearized bin width (double number)
 * peakPos        - Array of IRF peak positions (MxN double)
 * 
 * corrHist       - MxN uint32 matrix with resampled TDC histograms
 *
 * The calling syntax is:
 *
 *		corrHist = syntheticOJIP(dataStream, ...     % MxN uint16
 *                               firstCalBin, ...    % 1xN int16
 *                               lastCalBin, ...     % 1xN int16
 *                               realBinWidth, ...   % MxN double
 *                               linearBinWidth, ... % double
 *                               peakPos)            % MxN double
 *
 * Jakub Nedbal
 * King's College London
 * Aug 2021
 * Last Revision: 03-Aug-2021 - Make linearBinWidth a number not a vector
 *
 * Copyright 2021 Jakub Nedbal
 * BSD license
 *
 *=========================================================================
*/

#include "mex.h"
#include "math.h"
//#include "time.h"

#define	__dataStream        prhs[0]
#define	__firstCalBin       prhs[1]
#define	__lastCalBin    	prhs[2]
#define	__realBinWidth      prhs[3]
#define	__linearBinWidth    prhs[4]
#define	__peakPos           prhs[5]
#define	__outMatrix         plhs[0]

/* Validate inputs */
void validateInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 5 && nrhs != 6)
    { 
        mexErrMsgIdAndTxt("MATLAB:syntheticOJIP:invalidNumInputs",
                          "5 or 6 input arguments are required, ", 
                          "but %d were provided.", nrhs);
    }

    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:invalidNumOutputs",
            "1 output argument is required, but %d were asked for.", nlhs);
    }

    /* Ensure the inputs are of the correct data type */
    if (!mxIsUint32(__dataStream) && !mxIsInt32(__dataStream))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:invalidInputType",
            "1st input must be a 2D 32-bit integer matrix.");
    }

    if (!mxIsUint32(__firstCalBin) && !mxIsInt32(__firstCalBin))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:invalidInputType",
            "2nd input must be a 32-bit integer vector.");
    }

    if (!mxIsUint32(__lastCalBin) && !mxIsInt32(__lastCalBin))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:invalidInputType",
            "3rd input must be a 32-bit integer vector.");
    }

    if (!mxIsDouble(__realBinWidth) || mxIsComplex(__realBinWidth))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:invalidInputType",
            "4th input must be a 2D double matrix.");
    }

    if (!mxIsDouble(__linearBinWidth) || mxIsComplex(__linearBinWidth))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:invalidInputType",
            "5th input must be a double number.");
    }

    /* Check the array sizes */
    if (mxGetN(__dataStream) != mxGetNumberOfElements(__firstCalBin) || 
            mxGetN(__dataStream) != mxGetNumberOfElements(__lastCalBin) || 
            mxGetN(__dataStream) != mxGetN(__realBinWidth) || 
            mxGetNumberOfElements(__linearBinWidth) != 1)
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:wrongInputsSize",
            "2nd dimension of 1st input is not the same as "
            "number of elements in 2nd, 3rd input and "
            "the same as the 2nd dimension of the fourth input or"
            "the 5th input is not a single number.");
    }

//     /* Check the array sizes */
//     if (mxGetM(__dataStream) != mxGetM(__realBinWidth))
//     {
//         mexErrMsgIdAndTxt( "MATLAB:syntheticOJIP:wrongInputsSize",
//             "1st dimension of 1st input is not the same as "
//             "the 2nd dimension of the 4th input.");
//     }

    /* Specific check for optional 6th argument peak Pos */
    if (nrhs == 6)
    {
        if (!mxIsDouble(__peakPos) || mxIsComplex(__peakPos))
        {
            mexErrMsgIdAndTxt("MATLAB:syntheticOJIP:invalidInputType",
                "6th input must be a double vector.");
        }

        if (mxGetN(__dataStream) != mxGetNumberOfElements(__firstCalBin))
        {
            mexErrMsgIdAndTxt("MATLAB:syntheticOJIP:wrongInputsSize",
            "2nd dimension of 1st input is not the same as "
            "number of elements in 6th input and.");
        }
    }
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize numberPixels;    /* input scalar 
                       Number of pixels in the dataset */
    mwSize numberBins;      /* input scalar 
                       Number of TDC bins */
    mwSize numberFrames;    /* input scalar 
                       Number of data frames */
    int *outMatrix;    /* output matrix 
                       Matrix to hold the resampled histogram */
    int *dataStream;   /* input matrix 
                       Matrix of TDC histogram counts */
    int *firstCalBin;   /* input matrix 
                       Vector of first calibrated bins for each pixel */
    int *lastCalBin;    /* input matrix
                       Vector of last calibrated bins for each pixel */
    double *realBinWidth;   /* input matrix 
                       Matrix of actual TDC bin widths */
    double linearBinWidth;  /* input matrix 
                       Double number of linearized bin widths */
    double *peakPos;        /* input matrix 
                       Vector of corrected peak positions */
    double randMax;         /* number
                       Highest random number that can be produced */

    double randNr;          /* number
                       Number to store random code */
    double cumBinWidth;     /* number
                       Cummulative bin width for each pixel */
    double photonTime;      /* number
                       Simulated photon arrival time */
    int photonIndex;        /* number
                       Position of photon in linearized TDC histogram */
    int firstBin;           /* number
                       First calibrated bin in a pixel */
    int lastBin;            /* number
                       Last calibrated bin in a pixel */
    double pixelOffset;     /* number
                       Counting index of the first bin of each pixel */
    //clock_t t;              /* timer object */

    /* Validate the inputs */
    validateInputs(nlhs, plhs, nrhs, prhs);
    
   
    /* get the sizes of the input array */
    // Number of frames of the TDC
    numberFrames = mxGetM(__dataStream);
    mexPrintf("number of frames: %d\n", numberFrames);

    // Number of bins of the TDC
    numberBins = mxGetM(__realBinWidth);
    mexPrintf("number of bins: %d\n", numberBins);

    // Number of pixels in the sensor
    numberPixels = mxGetN(__dataStream);
    mexPrintf("number of pixels: %d\n", numberPixels);


    /* create pointers to the integer arrays in the input matrices  */
    dataStream = (int *) mxGetData(__dataStream);
    firstCalBin = (int *) mxGetData(__firstCalBin);
    lastCalBin = (int *) mxGetData(__lastCalBin);
    /* create pointers to the double arrays in the input matrices  */
    realBinWidth = mxGetPr(__realBinWidth);
    linearBinWidth = mxGetScalar(__linearBinWidth);
    /* Only if 6 input arguments are used */
    if (nrhs == 6)
    {
        peakPos = mxGetPr(__peakPos);
    }
    else
    {
        mxArray *__zeros;
        __zeros = mxCreateNumericMatrix((mwSize) numberPixels,
                                        1,
                                        mxDOUBLE_CLASS,
                                        mxREAL);
        peakPos = (double *) mxGetData(__zeros);
    }


    /* create the output matrix */
    __outMatrix = mxCreateNumericMatrix((mwSize) numberBins,
                                        (mwSize) numberFrames,
                                        mxINT32_CLASS,
                                        mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = (int *) mxGetData(__outMatrix);

    // Convert the maximum possible random number to a double
    randMax = (double) RAND_MAX;

    /* call the computational loop going pixel-by-pixel, bin-by-bin and 
       photon-by-photon */
    // Record the start time
    //t = clock();

//mexPrintf("ds: %d, %d, %d, %d\n", dataStream[0], dataStream[1], dataStream[2], dataStream[3]);
    // Go through a loop pixel by pixel
    for (int pixel = 0; pixel < numberPixels; pixel++)
    //for (int pixel = 0; pixel < 200; pixel++)
    {
        // Pixel offset
        //pixelOffset = pixel * numberBins;
        // Index of the first calibrated bin
        firstBin = firstCalBin[pixel];
        // Index of the last calibrated bin
        lastBin = lastCalBin[pixel];
        // Correct the pixel offset with the peak position
        pixelOffset = peakPos[pixel];
        //double time_taken = ((double) clock() - t) / CLOCKS_PER_SEC;
        mexPrintf("Pixel %d, elapsed time %f\n", pixel, 0);
        // Cumulative bin width
        cumBinWidth = 0;
        // Go through a loop frame by frame
        for (int bin = firstBin; bin < lastBin; bin++)
        //for (int bin = 3756; bin < 3757; bin++)
        {
            //mexPrintf("Pixel %d, bin %d, cumulative bin width %f\n", pixel, bin, cumBinWidth);
            //mexPrintf("ds: %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", dataStream[0], dataStream[1], dataStream[2], dataStream[3], dataStream[4], dataStream[5], dataStream[6], dataStream[7], dataStream[8], dataStream[9]);
            // Go through a loop frame by frame
            //for (int frame = 0; frame < numberFrames; frame++)
            for (int frame = 0; frame < numberFrames; frame++)
            {
                //mexPrintf("Pixel %d, bin %d, TDCvalue %d \n", frame * numberPixels + pixel, bin, dataStream[frame * numberPixels + pixel]);
                if (dataStream[frame * numberPixels + pixel] == bin)
                {
                    //mexPrintf("Index %d, Frame %d, bin %d\n", frame * numberPixels + pixel, frame, bin);
                    // Create a random number from 0 to 1
                    randNr = ((double) rand()) / randMax;
                    // mexPrintf("Random number %f, ", randNr);

                    // Scale the random number by the bin width and add it
                    // to the cummulative bin time
                    photonTime = realBinWidth[bin] * randNr + cumBinWidth;
                    // mexPrintf("Photon time %f, ", photonTime);

                    // Calculate the bin number of the simulated photon in
                    // the linearized histogram bins
                    photonIndex = (int) floor(photonTime / linearBinWidth + peakPos[pixel]);
                    // mexPrintf("Photon index %d, ", photonIndex);

                    // Add pixel offset to the pixel index
                    //photonIndex += pixelOffset;
                    // mexPrintf("Pixel corrected photon index %d\n", photonIndex);

                    // Update the output matrix
                    outMatrix[frame * numberBins + photonIndex] += 1;
                }
            }
            // Increment the cummulative bin width by the width of the
            // existing bin
            cumBinWidth += realBinWidth[bin];
        }
        // It is necessary to get rid of the last bin, as this will not be
        // complete and it would skew the data.
        // Find the index of the last bin
        photonIndex = (int) floor(cumBinWidth / linearBinWidth + peakPos[pixel]);
        // Add pixel offset to the bin index
        //photonIndex += pixelOffset;
        // Set this index to 0 to delete the conten of the last bin
        for (int frame = 0; frame < numberFrames; frame++)
        {
            outMatrix[photonIndex, frame] = 0;
        }    
        // mexPrintf("Pixel: %d\n", pixel);
    }
        
/*  mexPrintf("firstBin: %d\n", firstCalBin[0]);
    mexPrintf("lastBin: %d\n", lastCalBin[0]);
    mexPrintf("number of pixels: %d\n", numberPixels);
    mexPrintf("number of bins: %d\n", numberBins);
    mexPrintf("Value of firstBin: %d\n", inputHist[firstCalBin[0]]);
    mexPrintf("Value of lastBin: %d\n", inputHist[lastCalBin[0]]);
    mexPrintf("Maximum Rand Value: %d\n", RAND_MAX);*/
}