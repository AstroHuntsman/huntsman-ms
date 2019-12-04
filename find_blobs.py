#!/usr/bin/python3
from astropy.io import fits
import scipy.ndimage
import numpy as np
import sys
from scipy import stats

# Need to use special libraries to import and process FITS files

# Funky astro image processing
# import photutils.DAOStarFinder


# OpenCV

# Class to test whether to accept or reject blobs

class imframe():
    def __init__(self, diff, original, sigma):
        self.diff = diff
        self.original = original
        self.sigma = sigma

        # Flatten the arrays, and then create a corresponding meshgrid
        self.shape = self.diff.shape
        self.diff = self.diff.flatten()
        self.original = self.original.flatten()
        self.x = np.meshgrid(np.arange(diff.shape[0]), np.arange(diff.shape[1]), indexing='ij')
        self.x = [x.flatten() for x in self.x]
        self.flag = np.zeros(self.diff.size, dtype=np.int32)
        self.flag[np.nonzero(self.diff == 0)] = 1
        # added for blob_chi2_fast
        self.index = np.arange(self.diff.size)
        self.index.resize(self.shape)

    def blob_chi2(self, i, radius=5.0):
        x = self.x[0][i]
        y = self.x[1][i]
        xmin = np.int32(np.max([(x - radius), 0]))
        xmax = np.int32(np.min([(x + radius), self.shape[0]]))
        ymin = np.int32(np.max([(y - radius), 0]))
        ymax = np.int32(np.min([(y + radius), self.shape[1]]))
        inds = self.index[xmin:xmax, ymin:ymax].flatten().copy()
        r = np.sqrt(np.square(self.x[0][inds] - x) + np.square(self.x[1][inds] - y))
        indt = np.where(r < radius)[0]
        ind = inds[indt]
        r = r[indt]
        self.flag[ind] = 1

        # check the value of scale to use in stats.norm.pdf
        # Ideally radius=n*scale where n>3
        reduced_chi_squared = np.mean(np.square(self.diff[i] *
                                                stats.norm.pdf(r,
                                                               loc=0.0,
                                                               scale=1.0) - self.diff[ind]) /
                                      (self.sigma**2 + self.original[ind]))

        return(reduced_chi_squared)

    def findall(self, ind):
        self.lpeaks = []
        self.ipeaks = []
        indc = ind.copy()
        j = 0
        while indc.size > 0:
            print("j = {}".format(j))
            indc = indc[np.where(self.flag[indc] == 0)]
            if indc.size > 0:
                i = indc[np.argmax(self.diff[indc])]
                self.ipeaks.append(i)
                self.lpeaks.append(self.blob_chi2(i))
                j = j + 1
        for l, i in zip(self.lpeaks, self.ipeaks):
            print(f'chi^2: {l} at {i}')


def read_image(filename):
    fits.info(filename)

    # Now turn the image data into a numpy array
    raw_data = fits.getdata(filename)

    # Possibly want to convert to ints for speed
    # return(raw_data.astype(np.int32))
    return(raw_data.astype(np.float32))


def compare_data(raw1, raw2):

    # diff = ne.evaluate("raw2 - raw1")
    diff = raw2 - raw1

    hdu = fits.PrimaryHDU(data=diff)
    hdu.writeto("out_diff.fits", overwrite=True)
    return(diff)


def apply_laplacian(image_data):
    """ Calculates the Laplacian (grad-square) of a given array and saves it to file."""
    laplace = scipy.ndimage.laplace(image_data)

    # Save the data to a FITS file
    hdu = fits.PrimaryHDU(data=laplace)
    hdu.writeto("out_laplace.fits", overwrite=True)

    return(laplace)


def apply_gaussian(image_data):
    """ Calculates the Gaussian of a given array."""
    gaussian = scipy.ndimage.gaussian_filter(image_data, 4)

    # Save the data to a FITS file
    # hdu = fits.PrimaryHDU(data=gaussian)
    # hdu.writeto("out_gaussian.fits", overwrite=True)

    return(gaussian)


def get_pixels_above_threshold(data, thresh):
    """ Returns an array of the indices to all pixels with intensity greater than some threshold. """
    pixel_indices = np.nonzero(data > thresh)

    # Note: np.nonzero returns a two-element tuple of the row and column indices, so we need to take the
    # length of the first element in the tuple
    print("{} pixels greater than threshold".format(len(pixel_indices[0])))

    return(pixel_indices)


def main_from_files():

    # Get the filename from the script's arguments
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    raw1 = read_image(filename1)
    raw2 = read_image(filename2)

    diff = compare_data(raw1, raw2)

    # Define our threshold. The noise is artificially generated, so we can hard code the stddev
    stdev = 5.5105
    # thresh = 6*stdev
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


def main_from_stream():
    sys.path.append("/home/emily/Dropbox/Noodling/HuntsmanHack/huntsman-image-stream")
    from generate_images import ImageStream


def main_from_arrays(raw1, raw2):
    diff = compare_data(raw1, raw2)

    # Define our threshold. The noise is artificially generated, so we can hard code the stddev
    stdev = 5.5105
    # thresh = 6*stdev
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
    main_from_files()
