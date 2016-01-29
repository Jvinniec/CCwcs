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

#include <iterator>
#include "WcsHandler.h"

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
    char* cproj = str2char(proj) ;
    
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
    char* cctype1 = str2char(ctype1) ;
    char* cctype2 = str2char(ctype2) ;

    // Now do the wcs initialization
    wcs_ = wcskinit(naxis1, naxis2, cctype1, cctype2,
                    crpix1, crpix2, crval1, crval2, cd, cdelt1, cdelt2,
                    crota, equinox, epoch) ;
    
    delete[] cctype1 ;
    delete[] cctype2 ;
}


//_________________________________________________
// Create from a WorldCoor object
WcsHandler::WcsHandler(WorldCoor* wcs) :
    wcs_(wcs)
{}


//_________________________________________________
// Copy constructor
WcsHandler::WcsHandler(const WcsHandler& other) :
    wcs_(other.cloneWCS())
{}


//_________________________________________________
// Destructor
WcsHandler::~WcsHandler()
{
    // Free the memory that holds the wcs information
    wcsfree(wcs_) ;
}


# pragma mark - Public Methods

//_________________________________________________
// Create a new WorldCoor object with all the same values as wcs_
WorldCoor* WcsHandler::cloneWCS() const
{
    // Return a nullptr if this wcs_ is nullptr
    if (wcs_==nullptr) return nullptr ;
    
    // Otherwise we need to create an entirely new WorldCoor object
    WorldCoor* wcs_copy ;
    
    // Essentially all we want to do is to copy all of the values in the wcsobject
    wcs_copy->xref = wcs_->xref;		/* X reference coordinate value (deg) */
    wcs_copy->yref = wcs_->yref;		/* Y reference coordinate value (deg) */
    wcs_copy->xrefpix = wcs_->xrefpix;	/* X reference pixel */
    wcs_copy->yrefpix = wcs_->yrefpix;	/* Y reference pixel */
    wcs_copy->xinc = wcs_->xinc;		/* X coordinate increment (deg) */
    wcs_copy->yinc = wcs_->yinc;		/* Y coordinate increment (deg) */
    wcs_copy->rot = wcs_->rot;          /* rotation around axis (deg) (N through E) */
    /* rotation matrix */
    std::copy(std::begin(wcs_->cd), std::end(wcs_->cd), std::begin(wcs_copy->cd)) ;
    /* inverse rotation matrix */
    std::copy(std::begin(wcs_->dc), std::end(wcs_->dc), std::begin(wcs_copy->dc)) ;
    wcs_copy->equinox = wcs_->equinox;	/* Equinox of coordinates default to 1950.0 */
    wcs_copy->epoch = wcs_->epoch;		/* Epoch of coordinates default to equinox */
    wcs_copy->nxpix = wcs_->nxpix;		/* Number of pixels in X-dimension of image */
    wcs_copy->nypix = wcs_->nypix;		/* Number of pixels in Y-dimension of image */
    wcs_copy->plate_ra = wcs_->plate_ra;	/* Right ascension of plate center */
    wcs_copy->plate_dec = wcs_->plate_dec;	/* Declination of plate center */
    wcs_copy->plate_scale = wcs_->plate_scale;	/* Plate scale in arcsec/mm */
    wcs_copy->x_pixel_offset = wcs_->x_pixel_offset;	/* X pixel offset of image lower right */
    wcs_copy->y_pixel_offset = wcs_->y_pixel_offset;	/* Y pixel offset of image lower right */
    wcs_copy->x_pixel_size = wcs_->x_pixel_size;	/* X pixel_size */
    wcs_copy->y_pixel_size = wcs_->y_pixel_size;	/* Y pixel_size */
    /* pixel to plate coefficients for DSS */
    std::copy(std::begin(wcs_->ppo_coeff), std::end(wcs_->ppo_coeff), std::begin(wcs_copy->ppo_coeff)) ;
    /* X coefficients for plate model */
    std::copy(std::begin(wcs_->x_coeff), std::end(wcs_->x_coeff), std::begin(wcs_copy->x_coeff)) ;
    /* Y coefficients for plate model */
    std::copy(std::begin(wcs_->y_coeff), std::end(wcs_->y_coeff), std::begin(wcs_copy->y_coeff)) ;
    wcs_copy->xpix = wcs_->xpix;		/* X (RA) coordinate (pixels) */
    wcs_copy->ypix = wcs_->ypix;		/* Y (dec) coordinate (pixels) */
    wcs_copy->zpix = wcs_->zpix;		/* Z (face) coordinate (pixels) */
    wcs_copy->xpos = wcs_->xpos;		/* X (RA) coordinate (deg) */
    wcs_copy->ypos = wcs_->ypos;		/* Y (dec) coordinate (deg) */
    /* Values of CRPIXn keywords */
    std::copy(std::begin(wcs_->crpix), std::end(wcs_->crpix), std::begin(wcs_copy->crpix)) ;
    /* Values of CRVALn keywords */
    std::copy(std::begin(wcs_->crval), std::end(wcs_->crval), std::begin(wcs_copy->crval)) ;
    /* Values of CDELTn keywords */
    std::copy(std::begin(wcs_->cdelt), std::end(wcs_->cdelt), std::begin(wcs_copy->cdelt)) ;
    /* Values of PCiiijjj keywords */
    std::copy(std::begin(wcs_->pc), std::end(wcs_->pc), std::begin(wcs_copy->pc)) ;
    /* Constants for various projections */
    std::copy(std::begin(wcs_->projp), std::end(wcs_->projp), std::begin(wcs_copy->projp)) ;
    wcs_copy->pvfail = wcs_->pvfail;		/* If non-zero, significant inaccuracy likely to occur in projection */
    /* SCAMP constants for the PV coordinates */
    std::copy(std::begin(wcs_->projppv), std::end(wcs_->projppv), std::begin(wcs_copy->projppv)) ;
    
    wcs_copy->longpole = wcs_->longpole;	/* Longitude of North Pole in degrees */
    wcs_copy->latpole = wcs_->latpole;	/* Latitude of North Pole in degrees */
    wcs_copy->rodeg = wcs_->rodeg;		/* Radius of the projection generating sphere */
    wcs_copy->imrot = wcs_->imrot;		/* Rotation angle of north pole */
    wcs_copy->pa_north = wcs_->pa_north;	/* Position angle of north (0=horizontal) */
    wcs_copy->pa_east = wcs_->pa_east;	/* Position angle of east (0=horizontal) */
    wcs_copy->radvel = wcs_->radvel;		/* Radial velocity (km/sec away from observer)*/
    wcs_copy->zvel = wcs_->zvel;		/* Radial velocity (v/c away from observer)*/
    wcs_copy->zpzd = wcs_->zpzd;            /* Colat of FIP (degs) */
    wcs_copy->zpr = wcs_->zpr;              /* Radius of FIP (degs) */
    wcs_copy->imflip = wcs_->imflip;		/* If not 0, image is reflected around axis */
    wcs_copy->prjcode = wcs_->prjcode;      /* projection code (-1-32) */
    wcs_copy->latbase = wcs_->latbase;      /* Latitude base 90 (NPA), 0 (LAT), -90 (SPA) */
    wcs_copy->ncoeff1 = wcs_->ncoeff1;      /* Number of x-axis plate fit coefficients */
    wcs_copy->ncoeff2 = wcs_->ncoeff2;      /* Number of y-axis plate fit coefficients */
    wcs_copy->zpnp = wcs_->zpnp;            /* ZP polynomial order (0-9) */
    wcs_copy->changesys = wcs_->changesys;	/* 1 for FK4->FK5, 2 for FK5->FK4 */
                                            /* 3 for FK4->galactic, 4 for FK5->galactic */
    wcs_copy->printsys = wcs_->printsys;	/* 1 to print coordinate system, else 0 */
    wcs_copy->ndec = wcs_->ndec;            /* Number of decimal places in PIX2WCST */
    wcs_copy->degout = wcs_->degout;		/* 1 to always print degrees in PIX2WCST */
    wcs_copy->tabsys = wcs_->tabsys;		/* 1 to put tab between RA & Dec, else 0 */
    wcs_copy->rotmat = wcs_->rotmat;		/* 0 if CDELT, CROTA; 1 if CD */
    wcs_copy->coorflip = wcs_->coorflip;	/* 0 if x=RA, y=Dec; 1 if x=Dec, y=RA */
    wcs_copy->offscl = wcs_->offscl;		/* 0 if OK, 1 if offscale */
    wcs_copy->wcson = wcs_->wcson;          /* 1 if WCS is set, else 0 */
    wcs_copy->naxis = wcs_->naxis;          /* Number of axes in image (for WCSLIB 3.0) */
    wcs_copy->naxes = wcs_->naxes;		/* Number of axes in image */
    wcs_copy->wcsproj = wcs_->wcsproj;	/* WCS_OLD: AIPS worldpos() and worldpix()
                                           WCS_NEW: Mark Calabretta's WCSLIB subroutines
                                           WCS_BEST: WCSLIB for all but CAR,COE,NCP
                                           WCS_ALT:  AIPS for all but CAR,COE,NCP */
    wcs_copy->linmode = wcs_->linmode;	/* 0=system only, 1=units, 2=system+units */
    wcs_copy->detector = wcs_->detector;	/* Instrument detector number */
    /* Instrument name */
    std::copy(std::begin(wcs_->instrument), std::end(wcs_->instrument), std::begin(wcs_copy->instrument)) ;
    /* Values of CTYPEn keywords */
    std::copy(&wcs_->ctype[0][0], &wcs_->ctype[0][0]+(9*9), &wcs_copy->ctype[0][0]) ;
    /*  1st coordinate type code: RA--, GLON, ELON */
    std::copy(std::begin(wcs_->c1type), std::end(wcs_->c1type), std::begin(wcs_copy->c1type)) ;
    /*  2nd coordinate type code: DEC-, GLAT, ELAT */
    std::copy(std::begin(wcs_->c2type), std::end(wcs_->c2type), std::begin(wcs_copy->c2type)) ;
    /*  projection type code: SIN, TAN, ARC, NCP, GLS, MER, AIT, etc */
    std::copy(std::begin(wcs_->ptype), std::end(wcs_->ptype), std::begin(wcs_copy->ptype)) ;
    /* Units if LINEAR */
    std::copy(&wcs_->units[0][0], &wcs_->units[0][0]+(9*32), &wcs_copy->units[0][0]) ;
    /* Reference frame: FK4, FK4-NO-E, FK5, GAPPT*/
    std::copy(std::begin(wcs_->radecsys), std::end(wcs_->radecsys), std::begin(wcs_copy->radecsys)) ;
    /* Output reference frame: FK4,FK5,GAL,ECL */
    std::copy(std::begin(wcs_->radecout), std::end(wcs_->radecout), std::begin(wcs_copy->radecout)) ;
    /* Input reference frame: FK4,FK5,GAL,ECL */
    std::copy(std::begin(wcs_->radecin), std::end(wcs_->radecin), std::begin(wcs_copy->radecin)) ;
    wcs_copy->eqin = wcs_->eqin;		/* Input equinox (match sysin if 0.0) */
    wcs_copy->eqout = wcs_->eqout;		/* Output equinox (match sysout if 0.0) */
    wcs_copy->sysin = wcs_->sysin;		/* Input coordinate system code */
    wcs_copy->syswcs = wcs_->syswcs;	/* WCS coordinate system code */
    wcs_copy->sysout = wcs_->sysout;	/* Output coordinate system code */
                                        /* WCS_B1950, WCS_J2000, WCS_ICRS, WCS_GALACTIC,
                                         * WCS_ECLIPTIC, WCS_LINEAR, WCS_ALTAZ  */
    /* Center coordinates (with frame) */
    std::copy(std::begin(wcs_->center), std::end(wcs_->center), std::begin(wcs_copy->center)) ;
    wcs_copy->distcode = wcs_->distcode;	/* Distortion code 0=none 1=SIRTF */
    
    /* WCS command formats */
    /* where %s is replaced by WCS coordinates */
    /* where %f is replaced by the image filename */
    /* where %x is replaced by image coordinates */
    // Just copy them element by element
    for (int commandIndx=0; commandIndx<9; ++commandIndx)
        wcscominit(wcs_copy, commandIndx, wcs_->command_format[commandIndx]) ;
    
    /* Image rotation matrix */
    std::copy(std::begin(wcs_->ltm), std::end(wcs_->ltm), std::begin(wcs_copy->ltm)) ;
    /* Image offset */
    std::copy(std::begin(wcs_->ltv), std::end(wcs_->ltv), std::begin(wcs_copy->ltv)) ;
    /* First pixel to use in image (x, y) */
    std::copy(std::begin(wcs_->idpix), std::end(wcs_->idpix), std::begin(wcs_copy->idpix)) ;
    /* Number of pixels to use in image (x, y) */
    std::copy(std::begin(wcs_->ndpix), std::end(wcs_->ndpix), std::begin(wcs_copy->ndpix)) ;
    wcs_copy->wcsname = wcs_->wcsname;	/* WCS name (defaults to NULL pointer) */
    wcs_copy->wcschar = wcs_->wcschar;	/* WCS character (A-Z, null, space) */
    wcs_copy->logwcs = wcs_->logwcs;		/* 1 if DC-FLAG is set for log wavelength */
    
    // Handle structures
    wcs_copy->wcs = (WcsHandler(wcs_->wcs)).cloneWCS() ;
    wcs_copy->wcsdep = (WcsHandler(wcs_->wcsdep)).cloneWCS() ;
    
    // NOTE: I'll leave the following as an exercise to those who need them.
    //       If there's a reason that you do, create an issue in the gitHub repo
    //       and I'll work on getting it implemented.
//    struct poly	*inv_x;		/* SCAMP projection correction polynom in x */
//    struct poly	*inv_y;		/* SCAMP projection correction polynom in y */
//    struct Distort distort;	/* SIRTF distortion coefficients */
//    struct wcsprm wcsl;		/* WCSLIB main projection parameters */
//    struct linprm lin;		/* WCSLIB image/pixel conversion parameters */
//    struct celprm cel;		/* WCSLIB projection type */
//    struct prjprm prj;		/* WCSLIB projection parameters */
//    struct IRAFsurface *lngcor;	/* RA/longitude correction structure */
//    struct IRAFsurface *latcor;	/* Dec/latitude correction structure */
    
    // Return the fully copied object
    return wcs_copy ;
}


//_________________________________________________
/* wcsshift - Change the center of WCS */
void WcsHandler::ShiftImageCenter(double new_center_ra,     /* New center right ascension in degrees */
                                  double new_center_dec,	/* New center declination in degrees */
                                  const std::string& coordsys) /* FK4 or FK5 coordinates (1950 or 2000) */
{
    char* ccoordsys = str2char(coordsys) ;
    wcsshift(wcs_, new_center_ra, new_center_dec, ccoordsys) ;
    delete[] ccoordsys ;
}

//_________________________________________________
/* wcsfull - Return RA and Dec of image center, size in degrees */
void WcsHandler::GetRaDecSize_Deg(double  *cra,     /* Right ascension of image center (deg) (returned) */
                                  double  *cdec,	/* Declination of image center (deg) (returned) */
                                  double  *width,	/* Width in degrees (returned) */
                                  double  *height)  /* Height in degrees (returned) */
{
    wcsfull(wcs_, cra, cdec, width, height) ;
}

//_________________________________________________
/* wcsfull - Return RA and Dec of image center, size in radians */
void WcsHandler::GetRaDecSize_Rad(double  *cra,     /* Right ascension of image center (rad) (returned) */
                                  double  *cdec,	/* Declination of image center (rad) (returned) */
                                  double  *width,	/* Width in radians (returned) */
                                  double  *height)  /* Height in radians (returned) */
{
    // First get the values in degrees
    GetRaDecSize_Deg(cra, cdec, width, height) ;
    
    // Now convert all of the values to radians
    *cra *= D2R ;
    *cdec *= D2R ;
    *width *= D2R ;
    *height *= D2R ;
}


//_________________________________________________
/* wcssize - Return image center and size in RA and Dec */
void WcsHandler::GetImageCenterSize_Deg(double *cra,	/* Right ascension of image center (deg) (returned) */
                            double *cdec,	/* Declination of image center (deg) (returned) */
                            double *dra,	/* Half-width in right ascension (deg) (returned) */
                            double *ddec)	/* Half-width in declination (deg) (returned) */
{
    wcssize(wcs_, cra, cdec, dra, ddec) ;
}

//_________________________________________________
/* wcsrange - Return min and max RA and Dec of image in degrees */
void WcsHandler::GetImageRangeRaDec_Deg(double  *ra_min,	/* Min. right ascension of image (deg) (returned) */
                            double  *ra_max,	/* Max. right ascension of image (deg) (returned) */
                            double  *dec_min,	/* Min. declination of image (deg) (returned) */
                            double  *dec_max)	/* Max. declination of image (deg) (returned) */
{
    wcsrange(wcs_, ra_min, ra_max, dec_min, dec_max) ;
}


//_________________________________________________
/* wcscdset - Set scaling and rotation from CD matrix */
void WcsHandler::SetScalingRotation(std::vector<double> cd)
{
    if (cd.empty()) {
        wcscdset(wcs_, nullptr) ;
    } else if (cd.size() == 4) {
        wcscdset(wcs_, &cd[0]) ;
    } else {
        // Throw a runtime_error an error message
        throw std::runtime_error("[ERROR] WcsHandler::SetScalingRotation() : Provided vector must have size=4") ;
    }
}


//_________________________________________________
/* wcsdeltset - set scaling, rotation from CDELTi, CROTA2 */
void WcsHandler::SetScalingRotation(double cdelt1,	/* degrees/pixel in first axis (or both axes) */
                                    double cdelt2,	/* degrees/pixel in second axis if nonzero */
                                    double crota)	/* Rotation counterclockwise in degrees */
{
    wcsdeltset(wcs_, cdelt1, cdelt2, crota) ;
}


//_________________________________________________
/* wcspcset - set scaling, rotation from CDELTs and PC matrix */
void WcsHandler::SetScalingRotation(double cdelt1,	/* degrees/pixel in first axis (or both axes) */
                        double cdelt2,              /* degrees/pixel in second axis if nonzero */
                        std::vector<double> pc)     /* Rotation matrix, ignored if NULL */
{
    if (pc.empty()) {
        wcspcset(wcs_, cdelt1, cdelt2, nullptr) ;
    } else if (pc.size() == 4) {
        wcspcset(wcs_, cdelt1, cdelt2, &pc[0]) ;
    } else {
        // Throw a runtime_error an error message
        throw std::runtime_error("[ERROR] WcsHandler::SetScalingRotation() : Provided vector must have size=4") ;
    }
}


//_________________________________________________
/* setwcserr - Set WCS error message for later printing */
void WcsHandler::SetErrorMessage(const std::string& ErrMsg)
{
    char* cErrMsg = str2char(ErrMsg) ;
    setwcserr(cErrMsg) ;
    delete[] cErrMsg ;
}


//_________________________________________________
/* wcsoutinit - Set output coordinate system for pix2wcs */
void WcsHandler::SetWcsOutInit(const std::string& coord_sys)
/* coord_sys - Coordinate system (B1950, J2000, etc) */
{
    char* c_coord_sys = str2char(coord_sys) ;
    wcsoutinit(wcs_, c_coord_sys) ;
    delete[] c_coord_sys ;
}


//_________________________________________________
/* wcsininit - Set input coordinate system for wcs2pix */
void WcsHandler::WcsInInit(const std::string& coord_sys)
/* coord_sys - Coordinate system (B1950, J2000, etc) */
{
    char* c_coord_sys = str2char(coord_sys) ;
    wcsininit(wcs_, c_coord_sys) ;
    delete[] c_coord_sys ;
}


//_________________________________________________
/* wcsreset - Change WCS using arguments */
int WcsHandler::WcsReset(double crpix1,	/* Horizontal reference pixel */
             double crpix2,	/* Vertical reference pixel */
             double crval1,	/* Reference pixel horizontal coordinate in degrees */
             double crval2,	/* Reference pixel vertical coordinate in degrees */
             double cdelt1,	/* Horizontal scale in degrees/pixel, ignored if cd is not NULL */
             double cdelt2,	/* Vertical scale in degrees/pixel, ignored if cd is not NULL */
             double crota,	/* Rotation angle in degrees, ignored if cd is not NULL */
             std::vector<double> cd)	/* Rotation matrix, used if not NULL */
{
    if (cd.empty()) {
        return wcsreset(wcs_, crpix1, crpix2, crval1, crval2,
                        cdelt1, cdelt2, crota, nullptr) ;
    } else if (cd.size() == 4) {
        return wcsreset(wcs_, crpix1, crpix2, crval1, crval2,
                        cdelt1, cdelt2, crota, &cd[0]) ;
    } else {
        // Throw a runtime_error message
        throw std::runtime_error("[ERROR] WcsHandler::WcsReset() : Provided vector must have size=4") ;
        return -1 ;
    }
}


//_________________________________________________
/* wcscominit - Initialize catalog search command set by -wcscom */
void WcsHandler::WcsComInit(int i,                      /* Number of command (0-9) to initialize */
                            const std::string& command) /* command with %s where coordinates will go */
{
    char* com = str2char(command) ;
    wcscominit(wcs_, i, com) ;
    delete[] com ;
}


//_________________________________________________
/* wcscom - Execute catalog search command set by -wcscom */
void WcsHandler::WcsCom(int i,                          /* Number of command (0-9) to execute */
            const std::string& filename,    /* Image file name */
            double xfile,                   /* Horizontal image pixel coordinates for WCS command */
            double yfile,                   /* Vertical image pixel coordinates for WCS command */
            const std::string& wcstring)   /* WCS String from pix2wcst() */
{
    char* cfilename = str2char(filename) ;
    char* cwcstring = str2char(wcstring) ;
    wcscom(wcs_, i, cfilename, xfile, yfile, cwcstring) ;
    delete[] cfilename ;
    delete[] cwcstring ;
}



//_________________________________________________
/* savewcscom - Save WCS shell command */
void WcsHandler::SaveWcsCom(int i,                      /* i of 10 possible shell commands */
                            const std::string& wcscom)  /* Shell command using output WCS string */
{
    char* cwcscom = str2char(wcscom) ;
    savewcscom(i, cwcscom) ;
    delete[] cwcscom ;
}


//_________________________________________________
/* setwcsfile - Set filename for WCS error message */
void WcsHandler::SetWcsFile(const std::string& filename)
{
    char* cfilename = str2char(filename) ;
    setwcsfile(cfilename) ;
    delete[] cfilename ;
}



//_________________________________________________
/* cpwcs - Copy WCS keywords with no suffix to ones with suffix */
int WcsHandler::CpWcs(const std::string& header,    /* Pointer to start of FITS header */
                      const std::string& swcs)      /* Keyword suffix character for output WCS */
{
    char* cwcs = str2char(swcs) ;
    char* cheader = str2char(header) ;
    
    int ret = cpwcs(&cheader, cwcs) ;
    
    delete[] cwcs ;
    delete[] cheader ;
    
    return ret ;
}


//_________________________________________________
/* savewcscoor - Save output coordinate system */
void WcsHandler::SaveWcsCoor(const std::string& wcscoor) /* coordinate system (J2000, B1950, galactic) */
{
    char* cwcscoor = str2char(wcscoor) ;
    savewcscoor(cwcscoor) ;
    delete[] cwcscoor ;
}


//_________________________________________________
/* getwcscoor - Return output coordinate system */
std::string WcsHandler::GetWcsCoor()
{
    return std::string(getwcscoor()) ;
}


#pragma mark - Coordinate conversion subroutines in wcscon.c
/***********************************************
 * Coordinate conversion subroutines in wcscon.c
 ***********************************************/

//_________________________________________________
/* wcsconv - Convert between coordinate systems and equinoxes */
void WcsHandler::WcsConv(int coordsys_in,       /* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
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
                         double *rad_velocity)	/* Radial velocity in km/sec */
{
    // Access the underlying function
    wcsconv(coordsys_in, coordsys_out,
            equinox_in, equinox_out,
            epoch_in, epoch_out,
            x_coord, y_coord,
            x_coord_proper_motion, y_coord_proper_motion,
            parallax, rad_velocity) ;
}


//_________________________________________________
/* wcsconp - Convert between coordinate systems and equinoxes */
void WcsHandler::WcsConp(int coordsys_in,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
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
                         double *y_coord_proper_motion) /* Latitude or declination proper motion in degrees/year
                                                         Input in coordsys_in, returned in coordsys_out */
{
    // Access the underlying function (wcsconp)
    wcsconp(coordsys_in, coordsys_out,
            equinox_in, equinox_out,
            epoch_in, epoch_out,
            x_coord, y_coord,
            x_coord_proper_motion, y_coord_proper_motion) ;
}


//_________________________________________________
/* wcscon - Convert between coordinate systems and equinoxes */
void WcsHandler::WcsCon(int coordsys_in,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
                        int coordsys_out,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
                        double equinox_in,	/* Input equinox (default of coordsys_in if 0.0) */
                        double equinox_out, /* Output equinox (default of coordsys_out if 0.0) */
                        double *x_coord,	/* Longitude or right ascension in degrees
                                             Input in coordsys_in, returned in coordsys_out */
                        double *y_coord,	/* Latitude or declination in degrees
                                             Input in coordsys_in, returned in coordsys_out */
                        double epoch)       /* Besselian epoch in years */
{
    // Access the underlying function (wcscon)
    wcscon(coordsys_in, coordsys_out,
           equinox_in, equinox_out,
           x_coord, y_coord, epoch) ;
}


//_________________________________________________
/* fk425e - Convert B1950(FK4) to J2000(FK5) coordinates */
void WcsHandler::FK4to5e(double *ra,	/* Right ascension in degrees (B1950 in, J2000 out) */
                         double *dec,	/* Declination in degrees (B1950 in, J2000 out) */
                         double epoch)	/* Besselian epoch in years */
{
    // Access the underlying function (fk425e)
    fk425e(ra, dec, epoch) ;
}


//_________________________________________________
/* fk524e - Convert J2000(FK5) to B1950(FK4) coordinates */
void WcsHandler::FK5to4e(double *ra,	/* Right ascension in degrees (J2000 in, B1950 out) */
                         double *dec,	/* Declination in degrees (J2000 in, B1950 out) */
                         double epoch)	/* Besselian epoch in years */
{
    // Access the underlying function (fk524e)
    fk524e(ra, dec, epoch) ;
}


//_________________________________________________
/* wcscsys - Return code for coordinate system in string */
int WcsHandler::GetCoordSystemCode(const std::string& coord_sys) /* Coordinate system (B1950, J2000, etc) */
{
    char* ccoorsys = str2char(coord_sys) ;
    int code = wcscsys(ccoorsys) ;
    delete[] ccoorsys ;
    return code ;
}


//_________________________________________________
/* wcsceq - Set equinox from string (return 0.0 if not obvious) */
double WcsHandler::WcsCeq(const std::string& wcstring)  /* Coordinate system (B1950, J2000, etc) */
{
    char* cwcstring = str2char(wcstring) ;
    double ret = wcsceq(cwcstring);
    delete[] cwcstring ;
    return ret ;
}


//_________________________________________________
/* wcscstr - Set coordinate system type string from system and equinox */
std::string WcsHandler::WcsCstr(/* Coordinate system string (returned) */
                    int    coord_sys_code,	/* Coordinate system code */
                    double equinox,        /* Equinox of coordinate system */
                    double epoch)          /* Epoch of coordinate system */
{
    char* coord_sys ;
    wcscstr(coord_sys, coord_sys_code, equinox, epoch) ;
    return std::string(coord_sys) ;
}


//_________________________________________________
/* d2v3 - Convert RA and Dec in degrees and distance to vector */
/* x,y,z geocentric equatorial position of object (returned) */
std::vector<double> WcsHandler::D2v3(double	rra,	/* Right ascension in degrees */
                         double	rdec,	/* Declination in degrees */
                         double	r)      /* Distance to object in same units as pos */
{
    std::vector<double> pos(3,0.0) ;
    d2v3(rra, rdec, r, &pos[0]) ;
    return pos ;
}


//_________________________________________________
/* s2v3 - Convert RA and Dec in radians and distance to vector */
/* x,y,z geocentric equatorial position of object (returned) */
std::vector<double> WcsHandler::S2v3(double	rra,	/* Right ascension in radians */
                         double	rdec,	/* Declination in radians */
                         double	r)      /* Distance to object in same units as pos */
{
    std::vector<double> pos(3, 0.0) ;
    s2v3(rra, rdec, r, &pos[0]) ;
    return pos ;
}


//_________________________________________________
/* v2d3 - Convert vector to RA and Dec in degrees and distance */
void WcsHandler::V2d3(std::vector<double> pos,  /* x,y,z geocentric equatorial position of object */
          double	*rra,           /* Right ascension in degrees (returned) */
          double	*rdec,          /* Declination in degrees (returned) */
          double	*r)             /* Distance to object in same units as pos (returned) */
{
    v2d3(&pos[0], rra, rdec, r) ;
}


//_________________________________________________
/* v2s3 - Convert vector to RA and Dec in radians and distance */
void WcsHandler::V2s3(std::vector<double> pos,	/* x,y,z geocentric equatorial position of object */
                      double	*rra,           /* Right ascension in radians (returned) */
                      double	*rdec,          /* Declination in radians (returned) */
                      double	*r)             /* Distance to object in same units as pos (returned) */
{
    v2s3(&pos[0], rra, rdec, r) ;
}


#pragma mark - Distortion model subroutines in distort.c
/***********************************************
 * Distortion model subroutines in distort.c
 ***********************************************/

//_________________________________________________
/* distortinit - Set distortion coefficients from FITS header */
void WcsHandler::DistortInit(const std::string& hstring) /* FITS header */
{
    distortinit(wcs_, hstring.c_str()) ;
}


//_________________________________________________
/* setdistcode - Set WCS distortion code string from CTYPEi value */
void WcsHandler::SetDistCode(const std::string& header_ctype) /* CTYPE value from FITS header */
{
    char* c_header_ctype = str2char(header_ctype) ;
    setdistcode(wcs_, c_header_ctype) ;
    delete[] c_header_ctype ;
}


//_________________________________________________
/* DelDistort - Delete all distortion-related fields */
int WcsHandler::DeleteDistort(const std::string& header,    /* FITS header */
                  int verbose)                 /* If !=0, print keywords as deleted */
{
    char* cheader = str2char(header) ;
    int ret = DelDistort(cheader, verbose) ;
    delete[] cheader ;
    return ret ;
}


//_________________________________________________
/* pix2foc - Convert pixel to focal plane coordinates */
void WcsHandler::Pix2Foc(double x_pixel_coord,	/* Image pixel horizontal coordinate */
                         double y_pixel_coord,	/* Image pixel vertical coordinate */
                         double *x_focal_plane,	/* Focal plane horizontal coordinate(returned) */
                         double *y_focal_plane) /* Focal plane vertical coordinate (returned) */
{
    pix2foc(wcs_, x_pixel_coord, y_pixel_coord, x_focal_plane, y_focal_plane) ;
}


//_________________________________________________
/* foc2pix - Convert focal plane to pixel coordinates */
void WcsHandler::Foc2Pix(double x_focal_plane,	/* Focal plane horizontal coordinate(returned) */
                         double y_focal_plane,  /* Focal plane vertical coordinate (returned) */
                         double *x_pixel_coord,	/* Image pixel horizontal coordinate */
                         double *y_pixel_coord) /* Image pixel vertical coordinate */
{
    foc2pix(wcs_, x_focal_plane, y_focal_plane, x_pixel_coord, y_pixel_coord) ;
}


//_________________________________________________
/* SetFITSPlate - Set FITS header plate fit coefficients from structure */
/* Returns a string that represents the new header value */
std::string WcsHandler::SetFitsPlate(const std::string& header)  /* Current image FITS header */
{
    // First get a new header char* variable then fill it
    char* cheader = str2char(header) ;
    SetFITSPlate(cheader, wcs_) ;
    // Now create a new header string from the modified char* header and return it
    std::string new_header = std::string(cheader) ;
    delete[] cheader ;
    return new_header ;
}


//_________________________________________________
/* GetPlate - Return plate fit coefficients from structure in arguments */
int WcsHandler::GetPlateCoeff(int *ncoeff1,	/* Number of coefficients for x */
                              int *ncoeff2,	/* Number of coefficients for y) */
                              std::vector<double>& coeff) /* Plate fit coefficients */
{
    // I can cheat a little by knowing that the ncoeff1 and ncoeff2 come from wcs_
    // and that coeff is the size of the sum of the two
    int size = wcs_->ncoeff1 + wcs_->ncoeff2 ;
    
    // First create a coefficient array and fill it
    double* coeff_array = new double[size] ;
    int ret = GetPlate(wcs_, ncoeff1, ncoeff2, coeff_array) ;

    // Resize the coeff vector
    coeff.resize(*ncoeff1 + *ncoeff2) ;
    
    // Now fill the coeff vector from the values stored in the array
    for (int i=0; i < *ncoeff1; i++)
        coeff[i] = coeff_array[i] ;
    
    for (int i=0; i < *ncoeff2; i++)
        coeff[*ncoeff1+i] = coeff_array[*ncoeff1+1] ;
    
    delete[] coeff_array ;
    return ret ;
}

# pragma mark - Protected Methods

//_________________________________________________
/* Method for converting a string into a 'char*'. Note that this is a 'char*' which
 * is required by the methods in wcstools, not a 'const char*'. */
char* str2char(const std::string& str)
{
    char* cstr = new char[str.length()+1] ;
    std::strcpy(cstr, str.c_str()) ;
    return cstr ;
}
