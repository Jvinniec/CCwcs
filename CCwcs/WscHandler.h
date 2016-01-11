//
//  WscHandler.hpp
//  CCwcs
//
//  Created by Josh Cardenzana on 1/11/16.
//  Copyright © 2016 JCardenzana. All rights reserved.
//

#ifndef CCwcs_WscHandler_h
#define CCwcs_WscHandler_h

#include <stdio.h>
#include <string>
#include "wcs.h"

class WcsHandler {
public :
    // Default constructor
    WcsHandler() ;
    
    // Various constructors to reimplement the constructors in wcstools
    // set up WCS structure from a FITS image header
    WcsHandler(const std::string& header) ;
    // wcsxinit
    WcsHandler(double cra,      /* Center right ascension in degrees */
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
               const std::string& proj);	/* Projection */
    // wcskinit
    WcsHandler(int naxis1,	/* Number of pixels along x-axis */
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
               double epoch) ;
    
    // Copy constructor
    WcsHandler(const WcsHandler& other) ;
    
    // Destructor
    virtual ~WcsHandler() ;
    
    // Return a pointer to the contained WCS object
    WorldCoor* getWCS() {return wcs_ ;}
    
protected :
    // Keep track of the WCS information with an internal wcs object
    WorldCoor* wcs_ ;
    
private :
    
};

#endif /* WscHandler_h */
