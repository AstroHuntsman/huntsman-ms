#!/usr/bin/python3
import sys
sys.path.append("/home/emily/Dropbox/Noodling/HuntsmanHack/huntsman-image-stream")
from generate_images import ImageStream

import numpy as np
import numexpr as ne
import matplotlib.pyplot as plot
import scipy.ndimage
# Need to use special libraries to import and process FITS files
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

# Funky astro image processing
#import photutils.DAOStarFinder


# OpenCV
import cv2

# Class to test whether to accept or reject blobs
class imframe():
    def __init__(self,diff,original,sigma):
        self.diff=diff
        self.original = original
        self.sigma = sigma
    
        # Flatten the arrays, and then create a corresponding meshgrid
        self.shape=self.diff.shape
        self.diff=self.diff.flatten()
        self.original=self.original.flatten()
        self.x=np.meshgrid(np.arange(diff.shape[0]),np.arange(diff.shape[1]),indexing='ij')
        self.x=[x.flatten() for x in self.x]
        self.flag=np.zeros(self.diff.size,dtype=np.int32)
        self.flag[np.nonzero(self.diff == 0)] = 1

    def blob_chi2(self,i,radius=5.0):
        r=np.sqrt(np.square(self.x[0]-self.x[0][i])+np.square(self.x[1]-self.x[1][i]))
        # Pick some arbitrary (for now) radius of points to test
        ind=np.where(r<5)[0] 
        self.flag[ind]=1

        print("Indices in this test:\n{}".format(ind))
        print("Values at the indices:\n{}".format(self.diff[ind]))

        # TODO: Check this expression is correct
        return(np.mean(np.square(self.diff[i]*scipy.stats.norm.pdf(r[ind],loc=0.0,scale=5.0) -\
        self.diff[ind])/np.sqrt(self.sigma**2 + self.original[ind])))

    def findall(self,ind):
        self.lpeaks=[]
        self.ipeaks=[]
        indc=ind.copy()
        j=0        
        while indc.size > 0:
            print("j = {}".format(j))
            indc=indc[np.where(self.flag[indc]==0)]
            if indc.size > 0:
                i=indc[np.argmax(self.diff[indc])] 
                self.ipeaks.append(i)
                self.lpeaks.append(self.blob_chi2(i))
                j=j+1

def read_image(filename):
    fits.info(filename)

    # Now turn the image data into a numpy array
    raw_data = fits.getdata(filename)
    
    # Possibly want to convert to ints for speed
    #return(raw_data.astype(np.int32))
    return(raw_data.astype(np.float32))

def compare_data(raw1, raw2):

    #diff = ne.evaluate("raw2 - raw1")
    diff = raw2 - raw1

    hdu = fits.PrimaryHDU(data=diff)
    hdu.writeto("diff.fits", overwrite=True)
    return(diff)

def apply_laplacian(image_data):
    """ Calculates the Laplacian (grad-square) of a given array and saves it to file."""
    laplace = scipy.ndimage.laplace(image_data)

    # Save the data to a FITS file
    hdu = fits.PrimaryHDU(data=laplace)
    hdu.writeto("laplace.fits", overwrite=True)

    return(laplace)

def apply_gaussian(image_data):
    """ Calculates the Gaussian of a given array."""
    gaussian = scipy.ndimage.gaussian_filter(image_data, 4)

    # Save the data to a FITS file
    #hdu = fits.PrimaryHDU(data=gaussian)
    #hdu.writeto("gaussian.fits", overwrite=True)

    return(gaussian)

def get_pixels_above_threshold(data, thresh):
    """ Returns an array of the indices to all pixels with intensity greater than some threshold. """
    pixel_indices = np.nonzero(data > thresh)
    
    # Note: np.nonzero returns a two-element tuple of the row and column indices, so we need to take the
    # length of the first element in the tuple
    print("{} pixels greater than threshold".format(len(pixel_indices[0])))

    return(pixel_indices)

def main():

    # Get the filename from the script's arguments
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    raw1 = read_image(filename1)
    raw2 = read_image(filename2)
    
    diff = compare_data(raw1, raw2)

    # Define our threshold. The noise is artificially generated, so we can hard code the stddev
    stdev = 5.5105
    #thresh = 6*stdev
    # Set an arbitrary, but high, threshold for now. This is guaranteed to only accept pixels belonging
    # to the injected flash
    thresh = 1000

    print("Threshold = {}".format(thresh))
    data = np.resize(diff, diff.size)
    bright_indices = get_pixels_above_threshold(data, thresh)

    # Now calculate the likelihood of each blob being real
    # imframe needs the original data, as well as the differenced one. It also needs the standard
    # deviation
    im = imframe(diff, raw1, stdev)
    im.findall(bright_indices[0])
    print(im.lpeaks)

    # Finally, accept or reject the blobs, based on the chi-squared values
    # TODO

    
if __name__ == "__main__":
    main()