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
    
    /***********************************************
     * BELOW HERE ARE METHODS WHICH HAVE NOT CURRENTLY
     * BEEN IMPLEMENTED AS METHODS OF WcsHandler
     ***********************************************/
    
    void wcscom(	/* Execute catalog search command set by -wcscom */
                struct WorldCoor *wcs,	/* World coordinate system structure */
                int i,		/* Number of command (0-9) to execute */
                char *filename,	/* Image file name */
                double xfile,	/* Horizontal image pixel coordinates for WCS command */
                double yfile,	/* Vertical image pixel coordinates for WCS command */
                char *wcstring); /* WCS String from pix2wcst() */
    
    void savewcscom(	/* Save WCS shell command */
                    int i,		/* i of 10 possible shell commands */
                    char *wcscom);	/* Shell command using output WCS string */
    char *getwcscom(	/* Return WCS shell command */
                    int i);		/* i of 10 possible shell commands */
    void setwcscom(	/* Set WCS shell commands from stored values */
                   struct WorldCoor *wcs);	/* World coordinate system structure */
    void freewcscom(	/* Free memory storing WCS shell commands */
                    struct WorldCoor *wcs);	/* World coordinate system structure */
    
    void setwcsfile(	/* Set filename for WCS error message */
                    char *filename); /* FITS or IRAF file name */
    int cpwcs (		/* Copy WCS keywords with no suffix to ones with suffix */
               char **header,	/* Pointer to start of FITS header */
               char *cwcs);	/* Keyword suffix character for output WCS */
    
    void savewcscoor(	/* Save output coordinate system */
                     char *wcscoor);	/* coordinate system (J2000, B1950, galactic) */
    char *getwcscoor(void); /* Return output coordinate system */
    
    /* Coordinate conversion subroutines in wcscon.c */
    void wcsconv(	/* Convert between coordinate systems and equinoxes */
                 int sys1,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
                 int sys2,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
                 double eq1,	/* Input equinox (default of sys1 if 0.0) */
                 double eq2,	/* Output equinox (default of sys2 if 0.0) */
                 double ep1,	/* Input Besselian epoch in years */
                 double ep2,	/* Output Besselian epoch in years */
                 double *dtheta,	/* Longitude or right ascension in degrees
                                     Input in sys1, returned in sys2 */
                 double *dphi,	/* Latitude or declination in degrees
                                 Input in sys1, returned in sys2 */
                 double *ptheta,	/* Longitude or right ascension proper motion in deg/year
                                     Input in sys1, returned in sys2 */
                 double *pphi,	/* Latitude or declination proper motion in deg/year */
                 double *px,	/* Parallax in arcseconds */
                 double *rv);	/* Radial velocity in km/sec */
    void wcsconp(	/* Convert between coordinate systems and equinoxes */
                 int sys1,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
                 int sys2,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
                 double eq1,	/* Input equinox (default of sys1 if 0.0) */
                 double eq2,	/* Output equinox (default of sys2 if 0.0) */
                 double ep1,	/* Input Besselian epoch in years */
                 double ep2,	/* Output Besselian epoch in years */
                 double *dtheta,	/* Longitude or right ascension in degrees
                                     Input in sys1, returned in sys2 */
                 double *dphi,	/* Latitude or declination in degrees
                                 Input in sys1, returned in sys2 */
                 double *ptheta,	/* Longitude or right ascension proper motion in degrees/year
                                     Input in sys1, returned in sys2 */
                 double *pphi);	/* Latitude or declination proper motion in degrees/year
                                 Input in sys1, returned in sys2 */
    void wcscon(	/* Convert between coordinate systems and equinoxes */
                int sys1,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
                int sys2,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
                double eq1,	/* Input equinox (default of sys1 if 0.0) */
                double eq2,	/* Output equinox (default of sys2 if 0.0) */
                double *dtheta,	/* Longitude or right ascension in degrees
                                 Input in sys1, returned in sys2 */
                double *dphi,	/* Latitude or declination in degrees
                                 Input in sys1, returned in sys2 */
                double epoch);	/* Besselian epoch in years */
    void fk425e (	/* Convert B1950(FK4) to J2000(FK5) coordinates */
                 double *ra,	/* Right ascension in degrees (B1950 in, J2000 out) */
                 double *dec,	/* Declination in degrees (B1950 in, J2000 out) */
                 double epoch);	/* Besselian epoch in years */
    void fk524e (	/* Convert J2000(FK5) to B1950(FK4) coordinates */
                 double *ra,	/* Right ascension in degrees (J2000 in, B1950 out) */
                 double *dec,	/* Declination in degrees (J2000 in, B1950 out) */
                 double epoch);	/* Besselian epoch in years */
    int wcscsys(	/* Return code for coordinate system in string */
                char *coorsys);	 /* Coordinate system (B1950, J2000, etc) */
    double wcsceq (	/* Set equinox from string (return 0.0 if not obvious) */
                   char *wcstring);  /* Coordinate system (B1950, J2000, etc) */
    void wcscstr (	/* Set coordinate system type string from system and equinox */
                  char   *cstr,	 /* Coordinate system string (returned) */
                  int    syswcs,	/* Coordinate system code */
                  double equinox,	/* Equinox of coordinate system */
                  double epoch);	/* Epoch of coordinate system */
    void d2v3 (		/* Convert RA and Dec in degrees and distance to vector */
               double	rra,	/* Right ascension in degrees */
               double	rdec,	/* Declination in degrees */
               double	r,	/* Distance to object in same units as pos */
               double pos[3]);	/* x,y,z geocentric equatorial position of object (returned) */
    void s2v3 (		/* Convert RA and Dec in radians and distance to vector */
               double	rra,	/* Right ascension in radians */
               double	rdec,	/* Declination in radians */
               double	r,	/* Distance to object in same units as pos */
               double pos[3]);	/* x,y,z geocentric equatorial position of object (returned) */
    void v2d3 (		/* Convert vector to RA and Dec in degrees and distance */
               double	pos[3],	/* x,y,z geocentric equatorial position of object */
               double	*rra,	/* Right ascension in degrees (returned) */
               double	*rdec,	/* Declination in degrees (returned) */
               double	*r);	/* Distance to object in same units as pos (returned) */
    void v2s3 (		/* Convert vector to RA and Dec in radians and distance */
               double	pos[3],	/* x,y,z geocentric equatorial position of object */
               double	*rra,	/* Right ascension in radians (returned) */
               double	*rdec,	/* Declination in radians (returned) */
               double	*r);	/* Distance to object in same units as pos (returned) */
    
    /* Distortion model subroutines in distort.c */
    void distortinit (	/* Set distortion coefficients from FITS header */
                      struct WorldCoor *wcs,	/* World coordinate system structure */
                      const char* hstring);	/* FITS header */
    void setdistcode (	/* Set WCS distortion code string from CTYPEi value */
                      struct WorldCoor *wcs,	/* World coordinate system structure */
                      char	*ctype);	/* CTYPE value from FITS header */
    char *getdistcode (	/* Return distortion code string for CTYPEi */
                       struct WorldCoor *wcs);	/* World coordinate system structure */
    int DelDistort (	/* Delete all distortion-related fields */
                    char *header,	/* FITS header */
                    int verbose);	/* If !=0, print keywords as deleted */
    void pix2foc (	/* Convert pixel to focal plane coordinates */
                  struct WorldCoor *wcs,	/* World coordinate system structure */
                  double x,	/* Image pixel horizontal coordinate */
                  double y,	/* Image pixel vertical coordinate */
                  double *u,	/* Focal plane horizontal coordinate(returned) */
                  double *v);	/* Focal plane vertical coordinate (returned) */
    void foc2pix (	/* Convert focal plane to pixel coordinates */
                  struct WorldCoor *wcs,	/* World coordinate system structure */
                  double u,	/* Focal plane horizontal coordinate */
                  double v,	/* Focal plane vertical coordinate */
                  double *x,	/* Image pixel horizontal coordinate(returned) */
                  double *y);	/* Image pixel vertical coordinate (returned) */
    
    /* Other projection subroutines */
    
    /* 8 projections using AIPS algorithms (worldpos.c) */
    int worldpos (	/* Convert from pixel location to RA,Dec */
                  double xpix,	/* x pixel number  (RA or long without rotation) */
                  double ypix,	/* y pixel number  (Dec or lat without rotation) */
                  struct WorldCoor *wcs, /* WCS parameter structure */
                  double *xpos,	/* x (RA) coordinate (deg) (returned) */
                  double *ypos);	/* y (dec) coordinate (deg) (returned) */
    int worldpix (	/* Convert from RA,Dec to pixel location */
                  double xpos,	/* x (RA) coordinate (deg) */
                  double ypos,	/* y (dec) coordinate (deg) */
                  struct WorldCoor *wcs, /* WCS parameter structure */
                  double *xpix,	/* x pixel number (RA or long without rotation) */
                  double *ypix);	/* y pixel number (dec or lat without rotation) */
    
    /* Digital Sky Survey projection (dsspos.c) */
    int dsspos (	/* Convert from pixel location to RA,Dec */
                double xpix,	/* x pixel number  (RA or long without rotation) */
                double ypix,	/* y pixel number  (Dec or lat without rotation) */
                struct WorldCoor *wcs, /* WCS parameter structure */
                double *xpos,	/* x (RA) coordinate (deg) (returned) */
                double *ypos);	/* y (dec) coordinate (deg) (returned) */
    int dsspix (	/* Convert from RA,Dec to pixel location */
                double xpos,	/* x (RA) coordinate (deg) */
                double ypos,	/* y (dec) coordinate (deg) */
                struct WorldCoor *wcs, /* WCS parameter structure */
                double *xpix,	/* x pixel number (RA or long without rotation) */
                double *ypix);	/* y pixel number (dec or lat without rotation) */
    
    /* SAO TDC TAN projection with higher order terms (platepos.c) */
    int platepos (	/* Convert from pixel location to RA,Dec */
                  double xpix,	/* x pixel number  (RA or long without rotation) */
                  double ypix,	/* y pixel number  (Dec or lat without rotation) */
                  struct WorldCoor *wcs, /* WCS parameter structure */
                  double *xpos,	/* x (RA) coordinate (deg) (returned) */
                  double *ypos);	/* y (dec) coordinate (deg) (returned) */
    int platepix (	/* Convert from RA,Dec to pixel location */
                  double xpos,	/* x (RA) coordinate (deg) */
                  double ypos,	/* y (dec) coordinate (deg) */
                  struct WorldCoor *wcs, /* WCS parameter structure */
                  double *xpix,	/* x pixel number (RA or long without rotation) */
                  double *ypix);	/* y pixel number (dec or lat without rotation) */
    void SetFITSPlate (	/* Set FITS header plate fit coefficients from structure */
                       char *header,	/* Image FITS header */
                       struct WorldCoor *wcs); /* WCS structure */
    int SetPlate (	/* Set plate fit coefficients in structure from arguments */
                  struct WorldCoor *wcs, /* World coordinate system structure */
                  int ncoeff1,	/* Number of coefficients for x */
                  int ncoeff2,	/* Number of coefficients for y */
                  double *coeff);	/* Plate fit coefficients */
    int GetPlate (	/* Return plate fit coefficients from structure in arguments */
                  struct WorldCoor *wcs, /* World coordinate system structure */
                  int *ncoeff1,	/* Number of coefficients for x */
                  int *ncoeff2,	/* Number of coefficients for y) */
                  double *coeff);	/* Plate fit coefficients */
    
    /* IRAF TAN projection with higher order terms (tnxpos.c) */
    int tnxinit (	/* initialize the gnomonic forward or inverse transform */
                 const char *header, /* FITS header */
                 struct WorldCoor *wcs); /* pointer to WCS structure */
    int tnxpos (	/* forward transform (physical to world) gnomonic projection. */
                double xpix,	/* Image X coordinate */
                double ypix,	/* Image Y coordinate */
                struct WorldCoor *wcs, /* pointer to WCS descriptor */
                double *xpos,	/* Right ascension (returned) */
                double *ypos);	/* Declination (returned) */
    int tnxpix (	/* Inverse transform (world to physical) gnomonic projection */
                double xpos,     /* Right ascension */
                double ypos,     /* Declination */
                struct WorldCoor *wcs, /* Pointer to WCS descriptor */
                double *xpix,	/* Image X coordinate (returned) */
                double *ypix);	/* Image Y coordinate (returned) */
    
    /* IRAF ZPN projection with higher order terms (zpxpos.c) */
    int zpxinit (	/* initialize the zenithal forward or inverse transform */
                 const char *header, /* FITS header */
                 struct WorldCoor *wcs); /* pointer to WCS structure */
    int zpxpos (	/* forward transform (physical to world) */
                double xpix,	/* Image X coordinate */
                double ypix,	/* Image Y coordinate */
                struct WorldCoor *wcs, /* pointer to WCS descriptor */
                double *xpos,	/* Right ascension (returned) */
                double *ypos);	/* Declination (returned) */
    int zpxpix (	/* Inverse transform (world to physical) */
                double xpos,	/* Right ascension */
                double ypos,	/* Declination */
                struct WorldCoor *wcs, /* Pointer to WCS descriptor */
                double *xpix,	/* Image X coordinate (returned) */
                double *ypix);	/* Image Y coordinate (returned) */

    
protected :
    // Keep track of the WCS information with an internal wcs object
    WorldCoor* wcs_ ;
    
    // Method for converting a std::string into a char*
    char* str2char(const std::string& str) {
        char* cstr = new char[str.length()+1] ;
        return cstr ;
    }
    
    // Some internal helper methods based on the underlying values
    //double DegToRad() {return PI/180.0 ;}
    //double RadToDeg() {return 180.0/PI ;}
    
private :
    
};

#endif /* WscHandler_h */
