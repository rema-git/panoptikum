# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 23:03:24 2013

@author: Reinhardt A.W. Maier <rema@zaehlwerk.net>
"""
def singleImage(filename, part=[0,1392,0,1040]):
    """
    returns image data of a FITS primary HDU
    """
    import pyfits
    from numpy import rot90

    with pyfits.open(filename) as imagefile:
        imageContent = rot90(imagefile[0].data,3)[part[0]:part[1], part[2]:part[3]]
        return imageContent

def odTriple(filenames, part=[0,1392,0,1040], partLumCorr=[0,1392,0,1040], partLumCorrMask=[0,0,0,0]):
    """
    returns optical density calculated of three images (data of a FITS primary HDU)
    """
    from numpy import zeros, log

    mask = relative_to(partLumCorr, overlap(partLumCorr, partLumCorrMask))
    imageContent = zeros((3, part[1]-part[0], part[3]-part[2]))
    imageContentLumCorr = zeros(2)

    for n, filename in enumerate(filenames):
        imageContent[n] = singleImage(filename, part)
        if n < 2:
            LumCorrImage = singleImage(filename, partLumCorr)
            if mask is not None:
                LumCorrImage[mask[0]:mask[1], mask[2]:mask[3]] = 0
            imageContentLumCorr[n] = LumCorrImage.sum()

    # correct for wrong illumination of back image
    if (imageContentLumCorr != 0).all():
        imageContent[1] *= imageContentLumCorr[0] / imageContentLumCorr[1]

    # calculate optical density
    atomImage = imageContent[0] - imageContent[2]
    atomImage[atomImage<0] = 1

    backImage = imageContent[1] - imageContent[2]
    backImage[backImage<0] = 1

    od = -log(atomImage/backImage)

    return od

def overlap(r1, r2):
    """
    returns the overlapping part of two rectangles
    given by r1/r2 = [top, bottom, left, right]
    """
    if overlap_exists(r1, r2):
            olap = [sorted([max(r1[0],r2[0]), min(r1[1],r2[1])]), sorted([max(r1[2],r2[2]), min(r1[3],r2[3])])]
            return [y for x in olap for y in x]
    else:
        return False

def overlap_exists(r1, r2):
    """
    checks if two rectangles have an overlapping part
    """
    h_overlap = (r1[2] <= r2[3]) and (r1[3] >= r2[2])
    v_overlap = (r1[1] >= r2[0]) and (r1[0] <= r2[1])
    return h_overlap and v_overlap

def relative_to(r1, r2):
    """
    return r1 with respect to r2 coordinates
    """
    from numpy import size
    if size(r1) == 4 and size(r2) == 4:
        return [r2[0]-r1[0], r2[1]-r1[0], r2[2]-r1[2], r2[3]-r1[2]]

def WriteImage(image,filename):
    """
    write image to file in camera orientation
    """
    import pyfits
    from numpy import rot90

    image = rot90(image)
    hdu = pyfits.PrimaryHDU(image)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(filename, clobber=True)

def odTripleShowTarget(filenames, part, partLumCorr):
    """
    highlight in full frame image selected parts
    """
    import pylab
    pylab.clf()
    fullImage = odTriple(filenames)
    partImage = fullImage[part[0]:part[1], part[2]:part[3]]

    pylab.imshow(fullImage, vmin=partImage.min(), vmax=partImage.max())
    rect = pylab.Rectangle((part[2],part[0]), part[3]-part[2], part[1]-part[0],  ec='black', fc='none')
    rectLumCorr = pylab.Rectangle((partLumCorr[2],partLumCorr[0]), partLumCorr[3]-partLumCorr[2], partLumCorr[1]-partLumCorr[0],  ec='red', fc='none')
    pylab.gca().add_patch(rect)
    pylab.gca().add_patch(rectLumCorr)
