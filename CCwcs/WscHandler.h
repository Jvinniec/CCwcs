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
#include "wcs.h"

class WcsHandler {
public :
    // Default constructor
    WcsHandler() ;
    
    // Various constructors to reimplement the constructors in wcstools
    
    
    // Copy constructor
    
    // Destructor
    virtual ~WcsHandler() {} ;
    
    // Return a pointer to the contained WCS object
    WorldCoor* getWCS() {return wcs_ ;}
    
protected :
    // Keep track of the WCS information with an internal wcs object
    WorldCoor* wcs_ ;
    
private :
    
};

#endif /* WscHandler_h */
