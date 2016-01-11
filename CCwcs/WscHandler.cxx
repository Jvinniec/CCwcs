//
//  WscHandler.cxx
//  CCwcs
//
//  Created by Josh Cardenzana on 1/11/16.
//  Copyright © 2016 JCardenzana. All rights reserved.
//
//  Description: This class is meant to serve as a C++
//      wrapper for the wcstools which are written primarily
//      in C. In so doing the user doesnt need to worry about
//      deleting the objects within, or how to actually create
//      the object, only what a wcs object is meant to do.
//      Additionally, there are constructors derived from various
//      "wcsinit" methods to allow easy adoption with existing
//      code.

#include "WscHandler.h"

# pragma mark - Constructors/Destructors

//_________________________________________________
// Default constructor
WcsHandler::WcsHandler()
{}


//_________________________________________________
// Set up WCS structure from a FITS image header
WcsHandler::WcsHandler(const std::string& header) :
    wcs_(wcsinit(header.c_str()))
{}

//_________________________________________________
// wcsxinit
WcsHandler::WcsHandler(double cra,      /* Center right ascension in degrees */
                       double cdec,     /* Center declination in degrees */
                       double secpix,	/* Number of arcseconds per pixel */
                       double xrpix,	/* Reference pixel X coordinate */
                       double yrpix,	/* Reference pixel X coordinate */
                       int nxpix,       /* Number of pixels along x-axis */
                       int nypix,       /* Number of pixels along y-axis */
                       double rotate,	/* Rotation angle (clockwise positive) in degrees */
                       int equinox,     /* Equinox of coordinates, 1950 and 2000 supported */
                       double epoch,	/* Epoch of coordinates, used for FK4/FK5 conversion
                                         * no effect if 0 */
                       const std::string& proj) /* Projection */
{
    // Create the char* object
    char* cproj = new char[proj.length() + 1];
    std::strcpy(cproj, proj.c_str()) ;
    
    // Now do the wcs initialization
    wcs_ = wcsxinit(cra, cdec, secpix, xrpix, yrpix, nxpix, nypix,
                    rotate, equinox, epoch, cproj) ;
    
    delete[] cproj ;
}


//_________________________________________________
// wcskinit
WcsHandler::WcsHandler(int naxis1,	/* Number of pixels along x-axis */
                       int naxis2,	/* Number of pixels along y-axis */
                       const std::string& ctype1,	/* FITS WCS projection for axis 1 */
                       const std::string& ctype2,	/* FITS WCS projection for axis 2 */
                       double crpix1,	/* Reference pixel coordinates */
                       double crpix2,	/* Reference pixel coordinates */
                       double crval1,	/* Coordinate at reference pixel in degrees */
                       double crval2,	/* Coordinate at reference pixel in degrees */
                       double *cd,      /* Rotation matrix, used if not NULL */
                       double cdelt1,	/* scale in degrees/pixel, if cd is NULL */
                       double cdelt2,	/* scale in degrees/pixel, if cd is NULL */
                       double crota,	/* Rotation angle in degrees, if cd is NULL */
                       int equinox,     /* Equinox of coordinates, 1950 and 2000 supported */
                       double epoch) 	/* Epoch of coordinates, for FK4/FK5 conversion */
{
    // Create the char* objects
    char* cctype1 = new char[ctype1.length() + 1];
    char* cctype2 = new char[ctype2.length() + 1];
    std::strcpy(cctype1, ctype1.c_str()) ;
    std::strcpy(cctype2, ctype2.c_str()) ;
    
    // Now do the wcs initialization
    wcs_ = wcskinit(naxis1, naxis2, cctype1, cctype2,
                    crpix1, crpix2, crval1, crval2, cd, cdelt1, cdelt2,
                    crota, equinox, epoch) ;
    
    delete[] cctype1 ;
    delete[] cctype2 ;
}


//_________________________________________________
// Destructor
WcsHandler::~WcsHandler()
{
    // free the memory that holds the wcs information
    wcsfree(wcs_) ;
}


# pragma mark - Public Methods


# pragma mark - Protected Methods
