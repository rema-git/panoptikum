# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 20:40:53 2013

@author: Reinhardt A.W. Maier <rema@zaehlwerk.net>
"""
from numpy import size, argmax, exp, floor, sqrt, log, append, empty, nan, isnan, count_nonzero, pi, sort, abs
from scipy.optimize import curve_fit
from scipy import constants

def gauss(x, *p):
    A, mu, sigma, off = p
    return A * exp(-(x-mu)**2 / (2.*sigma**2)) + off

def linescan(image, dir=0, pixelSize=6.57e-6):
    return image.sum(axis=dir) * pixelSize

def linescanAtPos(image, position, dir=0):
    if dir == 1:
        return image[:,position]
    else:
        return image[position,:]

def gaussfit(image, dir=0, scale=1):
    try:
        linesum = linescan(image, dir)
        fit_initalVal = [max(linesum), argmax(linesum), 15, min(linesum)]
        coeff, var_matrix = curve_fit(gauss, range(0, size(image.conj().transpose(), dir)), linesum, p0=fit_initalVal)#, maxfev=50)
        if var_matrix.max() > 100:
            return empty(size(linesum)) * nan, empty(4) * nan, var_matrix
        fit = gauss(range(0,size(image.conj().transpose(), dir)), *coeff)
        return fit, coeff*scale, var_matrix
    except:
        return empty(size(linesum)) * nan , empty(4) * nan, empty((4, 4)) * nan

def atomCount(image, part=[0,0,0,0], pixelSize=6.57e-6, sigma=1.3565e-13):
    if part == [0,0,0,0]:
        part = [0,size(image,0),0,size(image,1)]
    try:
        atomCount = int(floor(pixelSize**2 / sigma * image[part[0]:part[1],part[2]:part[3]].sum()))
    except:
        atomCount = 0
    return atomCount

def atomCountGaussfit(gaussfitCoffH, gaussfitCoffV, pixelSize=6.57e-6, sigma=1.3565e-13):
    atomCount = empty(2)
    for m, n in enumerate([gaussfitCoffH, gaussfitCoffV]):
        Amp = n[1][0]
        width = n[1][2]
        atomCount[m] = Amp * sqrt(2*pi) * width * pixelSize / sigma
    if abs(atomCount[0] - atomCount[1]) > (atomCount[0] + atomCount[1]) / 4:
        print 'N_Gauss not well determined!'
        return nan
    return atomCount.mean()

def temperatur(tof, sigma, mass, insituSize=40e-6):
    if tof <= 0 or fwhm(sigma) < insituSize or fwhm(sigma) > 5e-3:
        return nan
    return mass*constants.m_u*(fwhm(sigma)**2 - insituSize**2) / (8*constants.k*log(2)*tof**2)

def fwhm(sigma, pixelSize=6.57e-6, reverse=0):
    if reverse:
        return abs(sigma / (2 * sqrt(2*log(2)) * pixelSize))
    else:
        return abs(2 * sqrt(2*log(2)) * sigma * pixelSize)

def analyseShot(shotPath, returnTemplate, image=[], filenames=['li_roi.fit', 'rb_roi.fit'], insituSizeLi=40e-6, insituSizeRb=40e-6):
    """
    returns measurment data of a single shot. Posible values for returnTemplate:
        * ID:            ADwin ID
        * ADwin varname: ADwin code value
        * N_Li/N_Rb:     atomnumber Li/Rb
        * e.g. fwhmH_Li: cloud width h/v
        * T_Li/T_Rb:     temperature
        * etc...
    """
    import os
    import adwin
    import apogee

    if image != []:
        image_Li = image
        image_Rb = image
        adwinFilename = shotPath
    else:
        if 'N_Li' in returnTemplate or 'fwhmH_Li' in returnTemplate or 'fwhmV_Li' in returnTemplate or 'T_Li' in returnTemplate:
            image_Li = apogee.singleImage(os.path.join(shotPath, filenames[0]))
        if 'N_Rb' in returnTemplate or 'fwhmH_Rb' in returnTemplate or 'fwhmV_Rb' in returnTemplate or 'T_Rb' in returnTemplate:
            image_Rb = apogee.singleImage(os.path.join(shotPath, filenames[1]))
        adwinFilename = os.path.join(shotPath, 'adwin_code.acm')

    if 'AmpH_Li' in returnTemplate or 'OffH_Li' in returnTemplate or 'PosH_Li' in returnTemplate or 'fwhmH_Li' in returnTemplate or 'N_Li_Gauss' in returnTemplate or 'T_Li' in returnTemplate:
       gaussfitCoffLiH = gaussfit(image_Li)
    if 'AmpV_Li' in returnTemplate or 'OffV_Li' in returnTemplate or 'PosV_Li' in returnTemplate or 'fwhmV_Li' in returnTemplate or 'N_Li_Gauss' in returnTemplate or 'T_Li' in returnTemplate:
       gaussfitCoffLiV = gaussfit(image_Li, dir = 1)
    if 'AmpH_Rb' in returnTemplate or 'OffH_Rb' in returnTemplate or 'PosH_Rb' in returnTemplate or 'fwhmH_Rb' in returnTemplate or 'N_Rb_Gauss' in returnTemplate or 'T_Rb' in returnTemplate:
       gaussfitCoffRbH = gaussfit(image_Rb)
    if 'AmpV_Rb' in returnTemplate or 'OffV_Rb' in returnTemplate or 'PosV_Rb' in returnTemplate or 'fwhmV_Rb' in returnTemplate or 'N_Rb_Gauss' in returnTemplate or 'T_Rb' in returnTemplate:
       gaussfitCoffRbV = gaussfit(image_Rb, dir = 1)

    returnList = empty([0, 0], dtype=object)

    for request in returnTemplate:
        if request == 'N_Li':
            returnList = append(returnList, atomCount(image_Li))
        elif request == 'N_Rb':
            returnList = append(returnList, atomCount(image_Rb))
        elif request == 'ID':
            returnList = append(returnList, os.path.split(shotPath)[-1])
        elif request == 'N_Rb_Gauss':
            returnList = append(returnList, atomCountGaussfit(gaussfitCoffRbH, gaussfitCoffRbV))
        elif request == 'N_Li_Gauss':
            returnList = append(returnList, atomCountGaussfit(gaussfitCoffLiH, gaussfitCoffLiV))
        elif request == 'fwhmH_Li':
            returnList = append(returnList, fwhm(gaussfitCoffLiH[1][2]))
        elif request == 'fwhmV_Li':
            returnList = append(returnList, fwhm(gaussfitCoffLiV[1][2]))
        elif request == 'fwhmH_Rb':
            returnList = append(returnList, fwhm(gaussfitCoffRbH[1][2]))
        elif request == 'fwhmV_Rb':
            returnList = append(returnList, fwhm(gaussfitCoffRbV[1][2]))
        elif request == 'OffH_Li':
            returnList = append(returnList, gaussfitCoffLiH[1][3])
        elif request == 'OffV_Li':
            returnList = append(returnList, gaussfitCoffLiV[1][3])
        elif request == 'OffH_Rb':
            returnList = append(returnList, gaussfitCoffRbH[1][3])
        elif request == 'OffV_Rb':
            returnList = append(returnList, gaussfitCoffRbV[1][3])
        elif request == 'PosH_Li':
            returnList = append(returnList, gaussfitCoffLiH[1][1])
        elif request == 'PosV_Li':
            returnList = append(returnList, gaussfitCoffLiV[1][1])
        elif request == 'PosH_Rb':
            returnList = append(returnList, gaussfitCoffRbH[1][1])
        elif request == 'PosV_Rb':
            returnList = append(returnList, gaussfitCoffRbV[1][1])
        elif request == 'AmpH_Li':
            returnList = append(returnList, gaussfitCoffLiH[1][0])
        elif request == 'AmpV_Li':
            returnList = append(returnList, gaussfitCoffLiV[1][0])
        elif request == 'AmpH_Rb':
            returnList = append(returnList, gaussfitCoffRbH[1][0])
        elif request == 'AmpV_Rb':
            returnList = append(returnList, gaussfitCoffRbV[1][0])
        elif request == 'T_Li':
            if abs(gaussfitCoffLiH[1][2] - gaussfitCoffLiV[1][2]) > (gaussfitCoffLiH[1][2]+gaussfitCoffLiV[1][2])/4:
                sigma = nan #max([gaussfitCoffLiV[1][2], gaussfitCoffLiH[1][2]])
                print 'T_Li not well determined!'
            else:
                sigma = (gaussfitCoffLiH[1][2]+gaussfitCoffLiV[1][2]) / 2
            returnList = append(returnList, temperatur(adwin.GetValue('ABB_TOF_LI', adwinFilename)*1e-6, sigma, 7, insituSizeLi))
        elif request == 'T_Rb':
            if abs(gaussfitCoffRbH[1][2] - gaussfitCoffRbV[1][2]) > (gaussfitCoffRbH[1][2]+gaussfitCoffRbV[1][2])/4:
                sigma = nan #max([gaussfitCoffRbV[1][2], gaussfitCoffRbH[1][2]])
                print 'T_Rb not well determined!'
            else:
                sigma = (gaussfitCoffRbH[1][2]+gaussfitCoffRbV[1][2]) / 2
            returnList = append(returnList, temperatur(adwin.GetValue('ABB_TOF_RB', adwinFilename)*1e-6, sigma, 87, insituSizeRb))
        else:
            returnList = append(returnList, adwin.GetValue(request, adwinFilename))

    return returnList
