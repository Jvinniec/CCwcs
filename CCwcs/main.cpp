//
//  main.cpp
//  CCwcs
//
//  Created by Josh Cardenzana on 1/11/16.
//  Copyright © 2016 JCardenzana. All rights reserved.
//

#include <iostream>

#include "WcsHandler.h"

int main(int argc, const char * argv[]) {
    // Note: This is only useful for the sake of testing
    
    /****************************************************
     * Attempt a conversion of celestial to galactic
     * coordinates as a test of this method.
     ****************************************************/
    
    // Setup some coordinates for tests, I like the Crab position
    double RA(83.6333), Dec(22.0145) ;
    double glon(RA), glat(Dec) ;
    
    std::cout << "Converting RA to galactic coordinates:" << std::endl;
    std::cout << "   RA : " << RA << std::endl;
    std::cout << "   Dec: " << Dec << std::endl;
    // Try to convert some basic coordinates from celestial to galactic coordinates
    WcsHandler::Celestial2Galactic(&glon, &glat) ;
    
    // Now print out the resulting computation
    std::cout << "Converted coordinates:" << std::endl;
    std::cout << "   GLon: " << glon << std::endl;
    std::cout << "   GLat: " << glat << std::endl;
    
    return 0 ;
}
