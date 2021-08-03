/*=========================================================================
 * syntheticPhotons.c - generating and binning synthetic photons
 *
 * This routine resamples histograms of photon arrival times
 * It takes the measured bin widths of actual TDC histograms.
 * It simulates randomly distributed photons within these histograms.
 * It then creates a new histogram with equally sized bins.
 *
 * Inputs: inputHist, firstCalBin, lastCalBin, realBinWidth, linearBinWidth
 * Output: corrHist
 *
 * inputHist      - Raw TDC histogram array (MxN uint32 array),
 *                  where M is number of bins, N is number of pixels
 * firstCalBin    - First indices of calibrated bins (1xN int32 vector)
 *                  Calibrated bins are those in the TDC range containing 
 *                  useful data
 * lastCalBin     - Last indices of calibrated bins (1xN int32 vector)
 * realBinWidth   - Array of calibrated bin widths (MxN double array)
 * linearBinWidth - Linearized bin width (double number)
 * peakPos        - Array of IRF peak positions (MxN double)
 * 
 * corrHist       - MxN uint32 matrix with resampled TDC histograms
 *
 * The calling syntax is:
 *
 *		corrHist = syntheticPhotons(inputHist, ...      % MxN uint32
 *                                  firstCalBin, ...    % 1xN int32
 *                                  lastCalBin, ...     % 1xN int32
 *                                  realBinWidth, ...   % MxN double
 *                                  linearBinWidth, ... % double
 *                                  peakPos)            % MxN double
 *
 * Jakub Nedbal
 * King's College London
 * May 2020
 * Last Revision: 15-Apr-2021 - Make linearBinWidth a number not a vector
 *
 * Copyright 2020 Jakub Nedbal
 * BSD license
 *
 *=========================================================================
*/

#include "mex.h"
#include "math.h"


#define	__inputHist         prhs[0]
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
        mexErrMsgIdAndTxt("MATLAB:syntheticPhotons:invalidNumInputs",
                          "5 or 6 input arguments are required, ", 
                          "but %d were provided.", nrhs);
    }

    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticPhotons:invalidNumOutputs",
            "1 output argument is required, but %d were asked for.", nlhs);
    }

    /* Ensure the inputs are of the correct data type */
    if (!mxIsUint32(__inputHist) && !mxIsInt32(__inputHist))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticPhotons:invalidInputType",
            "1st input must be a 2D 32-bit integer matrix.");
    }

    if (!mxIsUint32(__firstCalBin) && !mxIsInt32(__firstCalBin))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticPhotons:invalidInputType",
            "2nd input must be a 32-bit integer vector.");
    }

    if (!mxIsUint32(__lastCalBin) && !mxIsInt32(__lastCalBin))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticPhotons:invalidInputType",
            "3rd input must be a 32-bit integer vector.");
    }

    if (!mxIsDouble(__realBinWidth) || mxIsComplex(__realBinWidth))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticPhotons:invalidInputType",
            "4th input must be a 2D double matrix.");
    }

    if (!mxIsDouble(__linearBinWidth) || mxIsComplex(__linearBinWidth))
    {
        mexErrMsgIdAndTxt( "MATLAB:syntheticPhotons:invalidInputType",
            "5th input must be a double number.");
    }

    /* Check the array sizes */
    if (mxGetN(__inputHist) != mxGetNumberOfElements(__firstCalBin) || 
            mxGetN(__inputHist) != mxGetNumberOfElements(__lastCalBin) || 
            mxGetN(__inputHist) != mxGetN(__realBinWidth) || 
            mxGetNumberOfElements(__linearBinWidth) != 1)
    {
        mexErrMsgIdAndTxt( "MATLAB:max_in_place:wrongInputsSize",
            "2nd dimension of 1st input is not the same as "
            "number of elements in 2nd, 3rd input and "
            "the same as the 2nd dimension of the fourth input or"
            "the 5th input is not a single number.");
    }

    /* Check the array sizes */
    if (mxGetM(__inputHist) != mxGetM(__realBinWidth))
    {
        mexErrMsgIdAndTxt( "MATLAB:max_in_place:wrongInputsSize",
            "1st dimension of 1st input is not the same as "
            "the 2nd dimension of the 4th input.");
    }

    /* Specific check for optional 6th argument peak Pos */
    if (nrhs == 6)
    {
        if (!mxIsDouble(__peakPos) || mxIsComplex(__peakPos))
        {
            mexErrMsgIdAndTxt("MATLAB:syntheticPhotons:invalidInputType",
                "6th input must be a double vector.");
        }

        if (mxGetN(__inputHist) != mxGetNumberOfElements(__firstCalBin))
        {
            mexErrMsgIdAndTxt("MATLAB:max_in_place:wrongInputsSize",
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
    uint32_T *outMatrix;    /* output matrix 
                       Matrix to hold the resampled histogram */
    int *inputHist;         /* input matrix 
                       Matrix of TDC histogram counts */
    int *firstCalBin;       /* input matrix 
                       Vector of first calibrated bins for each pixel */
    int *lastCalBin;        /* input matrix
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
    int pixelOffset;        /* number
                       Counting index of the first bin of each pixel */

    /* Validate the inputs */
    validateInputs(nlhs, plhs, nrhs, prhs);
    
   
    /* get the sizes of the input array */
    // Number of bins of the TDC
    numberBins = mxGetM(__inputHist);

    // Number of pixels in the sensor
    numberPixels = mxGetN(__inputHist);
    

    /* create pointers to the integer arrays in the input matrices  */
    inputHist = (int *) mxGetData(__inputHist);
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
                                        (mwSize) numberPixels,
                                        mxUINT32_CLASS,
                                        mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = (uint32_T *) mxGetData(__outMatrix);
    
    // Convert the maximum possible random number to a double
    randMax = (double) RAND_MAX;

    /* call the computational loop going pixel-by-pixel, bin-by-bin and 
       photon-by-photon */

    // Go through a loop pixel by pixel
    for (int pixel = 0; pixel < numberPixels; pixel++)
    {
        // Pixel offset
        pixelOffset = pixel * numberBins;
        // Index of the first calibrated bin
        firstBin = pixelOffset + firstCalBin[pixel];
        // Index of the last calibrated bin
        lastBin = pixelOffset + lastCalBin[pixel];
        // Correct the pixel offset with the peak position
        pixelOffset += peakPos[pixel];
        // Cummulative bin width
        cumBinWidth = 0;
        // Go through a loop bin by bin
        for (int bin = firstBin; bin < lastBin; bin++)
        {
            // mexPrintf("Pixel %d, bin %d, cummulative bin width %f, number of photons %d\n", pixel, bin, cumBinWidth, inputHist[bin]);
            // Go through a loop photon by photon
            for (int photon = 0; photon < inputHist[bin]; photon++)
            {
                // Create a random number from 0 to 1
                randNr = ((double) rand()) / randMax;
                // mexPrintf("Random number %f, ", randNr);

                // Scale the random number by the bin width and add it to
                // the cummulative bin time
                photonTime = realBinWidth[bin] * randNr + cumBinWidth;
                // mexPrintf("Photon time %f, ", photonTime);

                // Calculate the bin number of the simulated photon in the
                // linearized histogram bins
                photonIndex = (int) floor(photonTime / linearBinWidth);
                // mexPrintf("Photon index %d, ", photonIndex);

                // Add pixel offset to the pixel index
                photonIndex += pixelOffset;
                // mexPrintf("Pixel corrected photon index %d\n", photonIndex);

                // Update the output matrix
                outMatrix[photonIndex] += 1;
            }
            // Increment the cummulative bin width by the width of the
            // existing bin
            cumBinWidth += realBinWidth[bin];
        }
        // It is necessary to get rid of the last bin, as this will
        // not be complete and it would skew the data.
        // Find the index of the last bin
        photonIndex = (int) floor(cumBinWidth / linearBinWidth);
        // Add pixel offset to the bin index
        photonIndex += pixelOffset;
        // Set this index to 0 to delete the conten of the last bin
        outMatrix[photonIndex] = 0;
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