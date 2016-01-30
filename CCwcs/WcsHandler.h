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
#include <vector>
#include "wcs.h"

class WcsHandler {
public :
    // Default constructor
    WcsHandler() ;
    
    // Various constructors to reimplement the constructors in wcstools.
    
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
               std::vector<double> cd,      /* Rotation matrix, used if not NULL */
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
    
    /**************************************************************************
     * Static methods for doing coordinate conversions
     **************************************************************************/
    
    // Some useful methods for converting coordinates between systems
    static void Celestial2Galactic(double* x_coord,  /* RA (deg) to be converted to Galactic Longitude (degrees) */
                                   double* y_coord,  /* Dec (deg) to be converted to Galactic Latitude (degrees) */
                                   double celestial_equinox=0.0,
                                   double galactic_equinox=0.0,
                                   double epoch=2000.0) ;
                                   
    static void Galactic2Celestial(double* x_coord,  /* Galactic Longitude (deg) to be converted to RA (degrees) */
                                   double* y_coord,  /* Galactic Latitude (deg) to be converted to Dec (degrees) */
                                   double celestial_equinox=0.0,
                                   double galactic_equinox=0.0,
                                   double epoch=2000.0) ;
    
    /* Note that in all Coordinate2Projected methods, the returned xproj,yprox are in degrees */
    static void Coordinate2Projected(double xcoord,                 /* Input Longitude/RA */
                                     double ycoord,                 /* Input Latitude/Dec */
                                     double xTangentPoint,          /* Longitude/RA Projection tangent point */
                                     double yTangentPoint,          /* Latitude/Dec Projection tangent point */
                                     double *xproj,                 /* Output x-projected coordinate */
                                     double *yproj) ;               /* Output y-projected coordinate */

    static void Coordinate2Projected(double xcoord,                 /* Input Longitude/RA */
                                     double ycoord,                 /* Input Latitude/Dec */
                                     double xTangentPoint,          /* Longitude/RA Projection tangent point */
                                     double yTangentPoint,          /* Latitude/Dec Projection tangent point */
                                     double *xproj,                 /* Output x-projected coordinate */
                                     double *yproj,                 /* Output y-projected coordinate */
                                     const std::string& coordSys,   /* Input coordinate system */
                                     const std::string& projection) ;   /* projection type */
    
    void Coordinate2Projected(double xcoord,    /* Input Longitude/RA */
                              double ycoord,    /* Input Latitude/Dec */
                              double *xproj,    /* Output x-projected coordinate */
                              double *yproj) ;  /* Output y-projected coordinate */
    
    /**************************************************************************
     * Static methods for calculating angular separation between two positions
     **************************************************************************/
    
    // Methods for getting the angular separation between two positions
    static double AngularSeparation_Deg(double x_coord1,
                                        double y_coord1,
                                        double x_coord2,
                                        double y_coord2);
    
    static double AngularSeparation_Rad(double x_coord1,
                                        double y_coord1,
                                        double x_coord2,
                                        double y_coord2) ;
    
    /**************************************************************************
     * Methods ported over from wcstools
     * Note that I've changed some of these names to make them a little more
     * meaningful. I've also identified the corresponding method in the wcstools
     **************************************************************************/
    
    /* wcstype - Set projection type from header CTYPEs */
    int SetWcsProjectionType(const std::string& ctype1,     /* FITS WCS projection for axis 1 */
                             const std::string& ctype2) ;   /* FITS WCS projection for axis 2 */
    
    /* iswcs - Returns 1 if wcs structure set, else 0 */
    int IsWcs() {return iswcs(wcs_) ;}
    
    /* nowcs - Returns 0 if wcs structure set, else 1 */
    int NoWcs() {return nowcs(wcs_) ;}
    
    /* pix2wcst - Convert pixel coordinates to World Coordinate string */
    int Pix2Wcst(double xpix,         /* Image horizontal coordinate in pixels */
                 double ypix,         /* Image vertical coordinate in pixels */
                 std::string& wcsstr, /* World coordinate string (returned) */
                 int lstr);           /* Length of world coordinate string (returned) */
    
    /* pix2wcs - Convert pixel coordinates to World Coordinates */
    void Pix2Wcs(double xpix,	/* Image horizontal coordinate in pixels */
                 double ypix,	/* Image vertical coordinate in pixels */
                 double *xpos,	/* Longitude/Right Ascension in degrees (returned) */
                 double *ypos)	/* Latitude/Declination in degrees (returned) */
        {return pix2wcs(wcs_, xpix, ypix, xpos, ypos) ;}
    
    /* wcsc2pix - Convert World Coordinates to pixel coordinates */
    void Wcsc2Pix(double xpos,	/* Longitude/Right Ascension in degrees */
                  double ypos,	/* Latitude/Declination in degrees */
                  const std::string& coorsys,/* Coordinate system (B1950, J2000, etc) */
                  double *xpix,	/* Image horizontal coordinate in pixels (returned) */
                  double *ypix,	/* Image vertical coordinate in pixels (returned) */
                  int *offscl);
    
    /* wcs2pix - Convert World Coordinates to pixel coordinates */
    void Wcs2Pix(double xpos,	/* Longitude/Right Ascension in degrees */
                 double ypos,	/* Latitude/Declination in degrees */
                 double *xpix,	/* Image horizontal coordinate in pixels (returned) */
                 double *ypix,	/* Image vertical coordinate in pixels (returned) */
                 int *offscl);
    
    
    /* wcsdist - Compute angular distance between 2 sky positions */
    static double WcsDist(double x_coord1,	/* First longitude/right ascension in degrees */
                          double y_coord1,	/* First latitude/declination in degrees */
                          double x_coord2,	/* Second longitude/right ascension in degrees */
                          double y_coord2)	/* Second latitude/declination in degrees */
        {return wcsdist(x_coord1, y_coord1, x_coord2, y_coord2) ;}
    
    /* wcsdist1 - Compute angular distance between 2 sky positions */
    static double WcsDist1(double x_coord1,	/* First longitude/right ascension in degrees */
                           double y_coord1,	/* First latitude/declination in degrees */
                           double x_coord2,	/* Second longitude/right ascension in degrees */
                           double y_coord2)	/* Second latitude/declination in degrees */
        {return wcsdist1(x_coord1, y_coord1, x_coord2, y_coord2) ;}
    
    /* wcsdiff - Compute angular distance between 2 sky positions */
    static double WcsDiff(double x_coord1,	/* First longitude/right ascension in degrees */
                          double y_coord1,	/* First latitude/declination in degrees */
                          double x_coord2,	/* Second longitude/right ascension in degrees */
                          double y_coord2)	/* Second latitude/declination in degrees */
        {return wcsdiff(x_coord1, y_coord1, x_coord2, y_coord2) ;}

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
    
    /* wcscdset - Set scaling and rotation from CD matrix */
    void SetScalingRotation(std::vector<double> cd=std::vector<double>(0)) ;
    
    /* wcsdeltset - set scaling, rotation from CDELTi, CROTA2 */
    void SetScalingRotation(double cdelt1,	/* degrees/pixel in first axis (or both axes) */
                            double cdelt2,	/* degrees/pixel in second axis if nonzero */
                            double crota);	/* Rotation counterclockwise in degrees */
    
    /* wcspcset - set scaling, rotation from CDELTs and PC matrix */
    void SetScalingRotation(double cdelt1,	/* degrees/pixel in first axis (or both axes) */
                            double cdelt2,	/* degrees/pixel in second axis if nonzero */
                            std::vector<double> pc=std::vector<double>(0));	/* Rotation matrix, ignored if NULL */
    
    /* setwcserr - Set WCS error message for later printing */
    void SetErrorMessage(const std::string& ErrMsg) ;
    
    /* wcserr - Print WCS error message to stderr */
    void PrintErrorMessage() {wcserr() ;}
    
    /* setdefwcs - Set flag to use AIPS WCS instead of WCSLIB */
    /* 1 for AIPS WCS subroutines, else WCSLIB */
    void SetAipsFlag(int oldwcs) {setdefwcs(oldwcs) ;}
    
    /* getdefwcs - Return flag for AIPS WCS set by setdefwcs */
    int GetDefWcs() {return getdefwcs() ;}
    
    /* getradecsys - Return name of image coordinate system */
    std::string GetRaDecSys() {return std::string(getradecsys(wcs_)) ;}
    
    /* wcsoutinit - Set output coordinate system for pix2wcs */
    void SetWcsOutInit(const std::string& coordsys) ;
    
    /* getwcsout - Return current output coordinate system */
    std::string GetWcsOut() {return std::string(getwcsout(wcs_) );}
    
    /* wcsininit - Set input coordinate system for wcs2pix */
    void WcsInInit(const std::string& coord_sys) ;
    
    /* getwcsin - Return current input coordinate system */
    std::string GetWcsIn() {return std::string(getwcsin(wcs_)) ;}
    
    /* setwcsdeg - Set WCS coordinate output format
     * degout : 1= degrees, 0= hh:mm:ss dd:mm:ss */
    int SetWcsDeg(int degout) { return setwcsdeg(wcs_, degout) ;}
    
    /* wcsndec - Set or get number of output decimal places
     * ndec : Number of decimal places in output string
     *        if < 0, return current ndec unchanged */
    int WcsNDec(int ndec) { return wcsndec(wcs_, ndec) ;}
    
    /* wcsreset - Change WCS using arguments */
    int WcsReset(double crpix1,	/* Horizontal reference pixel */
                 double crpix2,	/* Vertical reference pixel */
                 double crval1,	/* Reference pixel horizontal coordinate in degrees */
                 double crval2,	/* Reference pixel vertical coordinate in degrees */
                 double cdelt1,	/* Horizontal scale in degrees/pixel, ignored if cd is not NULL */
                 double cdelt2,	/* Vertical scale in degrees/pixel, ignored if cd is not NULL */
                 double crota,	/* Rotation angle in degrees, ignored if cd is not NULL */
                 std::vector<double> cd);	/* Rotation matrix, used if not NULL */
    
    /* wcseqset - Change equinox of reference pixel coordinates in WCS
     * equinox : Desired equinox as fractional year */
    void WcsEqSet(double equinox) {wcseqset(wcs_, equinox) ;}
    
    /* setwcslin - Set pix2wcst() mode for LINEAR coordinates
     * mode - 0: x y linear, 1: x units x units, 2: x y linear units */
    void SetWcsLin(int mode) { setwcslin(wcs_, mode) ;}
    
    /* wcszin - Set third dimension for cube projections
     * izpix - Set coordinate in third dimension (face) */
    int WcsZIn(int izpix) { return wcszin(izpix) ;}
    
    /* wcszout - Return coordinate in third dimension */
    int WcsZOut() { return wcszout(wcs_) ;}
    
    /* wcscominit - Initialize catalog search command set by -wcscom */
    void WcsComInit(int i,                          /* Number of command (0-9) to initialize */
                    const std::string& command);    /* command with %s where coordinates will go */
    
    /* wcscom - Execute catalog search command set by -wcscom */
    void WcsCom(int i,                          /* Number of command (0-9) to execute */
                const std::string& filename,    /* Image file name */
                double xfile,                   /* Horizontal image pixel coordinates for WCS command */
                double yfile,                   /* Vertical image pixel coordinates for WCS command */
                const std::string& wcstring) ;   /* WCS String from pix2wcst() */
    
    /* savewcscom - Save WCS shell command */
    void SaveWcsCom(int i,                       /* i of 10 possible shell commands */
                    const std::string& wcscom) ; /* Shell command using output WCS string */
    
    /* getwcscom - Return WCS shell command */
    std::string GetWcsCom(int i)                 /* i of 10 possible shell commands */
        {return std::string(getwcscom(i)) ;}
    
    /* setwcscom - Set WCS shell commands from stored values */
    void SetWcsCom() {setwcscom(wcs_) ;}
    
    /* freewcscom - Free memory storing WCS shell commands */
    void FreeWcsCom() {freewcscom(wcs_) ;}
    
    /* setwcsfile - Set filename for WCS error message */
    void SetWcsFile(const std::string& filename) ; /* FITS or IRAF file name */
    
    /* cpwcs - Copy WCS keywords with no suffix to ones with suffix */
    int CpWcs(const std::string& header,    /* Pointer to start of FITS header */
              const std::string& cwcs) ;    /* Keyword suffix character for output WCS */
    
    /* savewcscoor - Save output coordinate system */
    void SaveWcsCoor(const std::string& wcscoor) ; /* coordinate system (J2000, B1950, galactic) */

    /* getwcscoor - Return output coordinate system */
    std::string GetWcsCoor() ;
    
    /***********************************************
     * Coordinate conversion subroutines in wcscon.c
     ***********************************************/

    /* wcsconv - Convert between coordinate systems and equinoxes */
    static void WcsConv(int coordsys_in,       /* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
                 int coordsys_out,      /* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
                 double equinox_in,     /* Input equinox (default of coordsys_in if 0.0) */
                 double equinox_out,	/* Output equinox (default of coordsys_out if 0.0) */
                 double epoch_in,       /* Input Besselian epoch in years */
                 double epoch_out,      /* Output Besselian epoch in years */
                 double *x_coord,       /* Longitude or right ascension in degrees
                                         Input in sys1, returned in sys2 */
                 double *y_coord,       /* Latitude or declination in degrees
                                         Input in coordsys_in, returned in coordsys_out */
                 double *x_coord_proper_motion,/* Longitude or right ascension proper motion in deg/year
                                                Input in coordsys_in, returned in coordsys_out */
                 double *y_coord_proper_motion,	/* Latitude or declination proper motion in deg/year */
                 double *parallax,      /* Parallax in arcseconds */
                 double *rad_velocity);	/* Radial velocity in km/sec */
    
    /* wcsconp - Convert between coordinate systems and equinoxes */
    static void WcsConp(int coordsys_in,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
                 int coordsys_out,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
                 double equinox_in,	/* Input equinox (default of coordsys_in if 0.0) */
                 double equinox_out,/* Output equinox (default of coordsys_out if 0.0) */
                 double epoch_in,	/* Input Besselian epoch in years */
                 double epoch_out,	/* Output Besselian epoch in years */
                 double *x_coord,	/* Longitude or right ascension in degrees
                                     Input in coordsys_in, returned in coordsys_out */
                 double *y_coord,	/* Latitude or declination in degrees
                                     Input in coordsys_in, returned in coordsys_out */
                 double *x_coord_proper_motion,	/* Longitude or right ascension proper motion in degrees/year
                                                 Input in coordsys_in, returned in coordsys_out */
                 double *y_coord_proper_motion);/* Latitude or declination proper motion in degrees/year
                                                 Input in coordsys_in, returned in coordsys_out */
    
    /* wcscon - Convert between coordinate systems and equinoxes */
    static void WcsCon(int sys1,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
                int sys2,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
                double eq1,	/* Input equinox (default of sys1 if 0.0) */
                double eq2,	/* Output equinox (default of sys2 if 0.0) */
                double *dtheta,	/* Longitude or right ascension in degrees
                                 Input in sys1, returned in sys2 */
                double *dphi,	/* Latitude or declination in degrees
                                 Input in sys1, returned in sys2 */
                double epoch);	/* Besselian epoch in years */
    
    /* fk425e - Convert B1950(FK4) to J2000(FK5) coordinates */
    static void FK4to5e(double *ra,	/* Right ascension in degrees (B1950 in, J2000 out) */
                        double *dec,	/* Declination in degrees (B1950 in, J2000 out) */
                        double epoch);	/* Besselian epoch in years */
    
    /* fk524e - Convert J2000(FK5) to B1950(FK4) coordinates */
    static void FK5to4e(double *ra,	/* Right ascension in degrees (J2000 in, B1950 out) */
                        double *dec,	/* Declination in degrees (J2000 in, B1950 out) */
                        double epoch);	/* Besselian epoch in years */
    
    /* wcscsys - Return code for coordinate system in string */
    static int GetCoordSystemCode(const std::string& coord_sys) ; /* Coordinate system (B1950, J2000, etc) */
    
    /* wcsceq - Set equinox from string (return 0.0 if not obvious) */
    double WcsCeq(const std::string& wcstring) ; /* Coordinate system (B1950, J2000, etc) */
    
    /* wcscstr - Set coordinate system type string from system and equinox */
    std::string WcsCstr(/* Coordinate system string (returned) */
                        int    coord_sys_code,	/* Coordinate system code */
                        double equinox,        /* Equinox of coordinate system */
                        double epoch);         /* Epoch of coordinate system */
    
    /* d2v3 - Convert RA and Dec in degrees and distance to vector */
    /* x,y,z geocentric equatorial position of object (returned) */
    std::vector<double> D2v3(double	rra,	/* Right ascension in degrees */
                             double	rdec,	/* Declination in degrees */
                             double	r);     /* Distance to object in same units as pos */

    /* s2v3 - Convert RA and Dec in radians and distance to vector */
    /* x,y,z geocentric equatorial position of object (returned) */
    std::vector<double> S2v3(double	rra,	/* Right ascension in radians */
                             double	rdec,	/* Declination in radians */
                             double	r);     /* Distance to object in same units as pos */
    
    /* v2d3 - Convert vector to RA and Dec in degrees and distance */
    void V2d3(std::vector<double> pos,  /* x,y,z geocentric equatorial position of object */
              double	*rra,           /* Right ascension in degrees (returned) */
              double	*rdec,          /* Declination in degrees (returned) */
              double	*r);            /* Distance to object in same units as pos (returned) */
    
    /* v2s3 - Convert vector to RA and Dec in radians and distance */
    void V2s3(std::vector<double> pos,	/* x,y,z geocentric equatorial position of object */
              double	*rra,           /* Right ascension in radians (returned) */
              double	*rdec,          /* Declination in radians (returned) */
              double	*r);            /* Distance to object in same units as pos (returned) */
    
    /***********************************************
     * Distortion model subroutines in distort.c
     ***********************************************/

    /* distortinit - Set distortion coefficients from FITS header */
    void DistortInit(const std::string& hstring); /* FITS header */

    /* setdistcode - Set WCS distortion code string from CTYPEi value */
    void SetDistCode(const std::string& header_ctype) ; /* CTYPE value from FITS header */

    /* getdistcode - Return distortion code string for CTYPEi */
    std::string GetDistCode() {return std::string(getdistcode(wcs_)) ;}
    
    /* DelDistort - Delete all distortion-related fields */
    int DeleteDistort(const std::string& header,    /* FITS header */
                      int verbose) ;                /* If !=0, print keywords as deleted */
    
    /* pix2foc - Convert pixel to focal plane coordinates */
    void Pix2Foc(double x_pixel_coord,	/* Image pixel horizontal coordinate */
                 double y_pixel_coord,	/* Image pixel vertical coordinate */
                 double *x_focal_plane,	/* Focal plane horizontal coordinate(returned) */
                 double *y_focal_plane);/* Focal plane vertical coordinate (returned) */
    
    /* foc2pix - Convert focal plane to pixel coordinates */
    void Foc2Pix(double x_focal_plane,	/* Focal plane horizontal coordinate(returned) */
                 double y_focal_plane,  /* Focal plane vertical coordinate (returned) */
                 double *x_pixel_coord,	/* Image pixel horizontal coordinate */
                 double *y_pixel_coord);/* Image pixel vertical coordinate */
    
    /***********************************************
     * Other projection subroutines
     ***********************************************/
    
    /* 8 projections using AIPS algorithms (worldpos.c) */
    
    /* worldpos - Convert from pixel location to RA,Dec */
    int WorldPos(double xpix,	/* x pixel number  (RA or long without rotation) */
                 double ypix,	/* y pixel number  (Dec or lat without rotation) */
                 double *xpos,	/* x (RA) coordinate (deg) (returned) */
                 double *ypos)	/* y (dec) coordinate (deg) (returned) */
        {return worldpos(xpix, ypix, wcs_, xpos, ypos) ;}
    
    /* worldpix - Convert from RA,Dec to pixel location */
    int WorldPix(double xpos,	/* x (RA) coordinate (deg) */
                 double ypos,	/* y (dec) coordinate (deg) */
                 double *xpix,	/* x pixel number (RA or long without rotation) */
                 double *ypix)	/* y pixel number (dec or lat without rotation) */
        {return worldpix(xpos, ypos, wcs_, xpix, ypix) ;}
    
    /***********************************************
     * Digital Sky Survey projection (dsspos.c)
     ***********************************************/
    
    /* dsspos - Convert from pixel location to RA,Dec */
    int DssPos(double xpix,     /* x pixel number  (RA or long without rotation) */
               double ypix,     /* y pixel number  (Dec or lat without rotation) */
               double *xpos,	/* x (RA) coordinate (deg) (returned) */
               double *ypos)	/* y (dec) coordinate (deg) (returned) */
        {return dsspos(xpix, ypix, wcs_, xpos, ypos) ;}
    
    /* dsspix - Convert from RA,Dec to pixel location */
    int DssPix(double xpos,     /* x (RA) coordinate (deg) */
               double ypos,     /* y (dec) coordinate (deg) */
               double *xpix,	/* x pixel number (RA or long without rotation) */
               double *ypix)	/* y pixel number (dec or lat without rotation) */
        {return dsspix(xpos, ypos, wcs_, xpix, ypix) ;}
    
    /***********************************************
     * SAO TDC TAN projection with higher order terms (platepos.c)
     ***********************************************/
    
    /* platepos - Convert from pixel location to RA,Dec */
    int PlatePos(double xpix,	/* x pixel number  (RA or long without rotation) */
                 double ypix,	/* y pixel number  (Dec or lat without rotation) */
                 double *xpos,	/* x (RA) coordinate (deg) (returned) */
                 double *ypos)	/* y (dec) coordinate (deg) (returned) */
        {return platepos(xpix, ypix, wcs_, xpos, ypos) ;}
    
    /* platepix - Convert from RA,Dec to pixel location */
    int PlatePix(double xpos,	/* x (RA) coordinate (deg) */
                 double ypos,	/* y (dec) coordinate (deg) */
                 double *xpix,	/* x pixel number (RA or long without rotation) */
                 double *ypix)	/* y pixel number (dec or lat without rotation) */
        {return platepix(xpos, ypos, wcs_, xpix, ypix) ;}

    /* SetFITSPlate - Set FITS header plate fit coefficients from structure */
    /* Returns a string that represents the new header value */
    std::string SetFitsPlate(const std::string& header) ; /* Current image FITS header */
    
    /* SetPlate - Set plate fit coefficients in structure from arguments */
    int SetPlateCoeff(int ncoeff1,	/* Number of coefficients for x */
                      int ncoeff2,	/* Number of coefficients for y */
                      std::vector<double> coeff)	/* Plate fit coefficients */
        {return SetPlate(wcs_, ncoeff1, ncoeff2, &coeff[0]) ;}
    
    /* GetPlate - Return plate fit coefficients from structure in arguments */
    int GetPlateCoeff(int *ncoeff1,	/* Number of coefficients for x */
                      int *ncoeff2,	/* Number of coefficients for y) */
                      std::vector<double>& coeff) ; /* Plate fit coefficients */
    
    /***********************************************
     * IRAF TAN projection with higher order terms (tnxpos.c)
     ***********************************************/
    
    /* tnxinit - initialize the gnomonic forward or inverse transform */
    int TnxInit(const std::string& header) /* FITS header */
        {return tnxinit(header.c_str(), wcs_) ;}
    
    /* tnxpos - forward transform (physical to world) gnomonic projection. */
    int TnxPos(double xpix,     /* Image X coordinate */
               double ypix,     /* Image Y coordinate */
               double *xpos,	/* Right ascension (returned) */
               double *ypos)	/* Declination (returned) */
        {return tnxpos(xpix, ypix, wcs_, xpos, ypos) ;}
    
    /* tnxpix - Inverse transform (world to physical) gnomonic projection */
    int TnxPix(double xpos,     /* Right ascension */
               double ypos,     /* Declination */
               double *xpix,	/* Image X coordinate (returned) */
               double *ypix)	/* Image Y coordinate (returned) */
        {return tnxpix(xpos, ypos, wcs_, xpix, ypix) ;}
    
    /***********************************************
     * IRAF ZPN projection with higher order terms (zpxpos.c)
     ***********************************************/
    
    /* zpxinit - initialize the zenithal forward or inverse transform */
    int ZpxInit(const std::string& header)
        {return zpxinit(header.c_str(), wcs_) ;}
    
    /* zpxpos - forward transform (physical to world) */
    int ZpxPos(double xpix,     /* Image X coordinate */
               double ypix,     /* Image Y coordinate */
               double *xpos,	/* Right ascension (returned) */
               double *ypos)	/* Declination (returned) */
        {return zpxpos(xpix, ypix, wcs_, xpos, ypos) ;}
    
    /* zpxpix - Inverse transform (world to physical) */
    int ZpxPix(double xpos,     /* Right ascension */
               double ypos,     /* Declination */
               double *xpix,	/* Image X coordinate (returned) */
               double *ypix)	/* Image Y coordinate (returned) */
        {return zpxpix(xpos, ypos, wcs_, xpix, ypix) ;}
    
protected :
    // Keep track of the WCS information with an internal wcs object
    WorldCoor* wcs_ ;
    
    // Method for converting a std::string into a char*
    static char* str2char(const std::string& str) ;
    
private :
    
};

#endif /* WscHandler_h */
