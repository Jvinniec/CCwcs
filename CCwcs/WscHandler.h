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
    
    // Create from a WorldCoor object
    WcsHandler(WorldCoor* wcs) ;
    // Copy constructor
    WcsHandler(const WcsHandler& other) ;
    
    // Destructor
    virtual ~WcsHandler() ;
    
    // Return a pointer to the contained WCS object
    WorldCoor* getWCS() {return wcs_ ;}
    // Return a separate copy of the contained WCS object
    WorldCoor* cloneWCS() const ;
    
    // Methods ported over from wcstools ===============================
    // Note that I've changed some of these names to make them a little more
    // meaningful. I've also identified the corresponding method in the wcstools
    
    /* wcsshift - Change center of WCS */
    void ShiftImageCenter(double new_center_ra,     /* New center right ascension in degrees */
                          double new_center_dec,	/* New center declination in degrees */
                          const std::string& coordsys); /* FK4 or FK5 coordinates (1950 or 2000) */
    
    /* wcsfull - Return RA and Dec of image center, size in degrees */
    void GetRaDecSize_Deg(double  *cra,     /* Right ascension of image center (deg) (returned) */
                          double  *cdec,	/* Declination of image center (deg) (returned) */
                          double  *width,	/* Width in degrees (returned) */
                          double  *height) ;/* Height in degrees (returned) */
    void GetRaDecSize_Rad(double  *cra,     /* Right ascension of image center (rad) (returned) */
                          double  *cdec,	/* Declination of image center (rad) (returned) */
                          double  *width,	/* Width in radians (returned) */
                          double  *height) ;/* Height in radians (returned) */
    
    /* wcscent - Print the image center and size in WCS units */
    void PrintCenter() {wcscent(wcs_) ;}
    
    /* wcssize - Return image center and size in RA and Dec */
    void GetImageCenterSize_Deg(double *cra,	/* Right ascension of image center (deg) (returned) */
                                double *cdec,	/* Declination of image center (deg) (returned) */
                                double *dra,	/* Half-width in right ascension (deg) (returned) */
                                double *ddec);	/* Half-width in declination (deg) (returned) */
    
    /* wcsrange - Return min and max RA and Dec of image in degrees */
    void GetImageRangeRaDec_Deg(double  *ra_min,	/* Min. right ascension of image (deg) (returned) */
                                double  *ra_max,	/* Max. right ascension of image (deg) (returned) */
                                double  *dec_min,	/* Min. declination of image (deg) (returned) */
                                double  *dec_max);	/* Max. declination of image (deg) (returned) */
    
protected :
    // Keep track of the WCS information with an internal wcs object
    WorldCoor* wcs_ ;
    
    // Some internal helper methods based on the underlying values
    //double DegToRad() {return PI/180.0 ;}
    //double RadToDeg() {return 180.0/PI ;}
    
private :
    
};

#endif /* WscHandler_h */
