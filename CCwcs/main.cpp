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
    
    std::cout << "Converting celestial to Galactic coordinates:" << std::endl;
    std::cout << "   RA : " << RA << " (degrees)" << std::endl;
    std::cout << "   Dec: " << Dec << " (degrees)" << std::endl;
    // Try to convert some basic coordinates from celestial to galactic coordinates
    WcsHandler::Celestial2Galactic(&glon, &glat) ;
    
    // Now print out the resulting computation
    std::cout << "   Converted coordinates:" << std::endl;
    std::cout << "   GLon: " << glon << " (degrees)" << std::endl;
    std::cout << "   GLat: " << glat << " (degrees)" << std::endl;
    
    /****************************************************
     * Attempt a conversion of galactic to celestial
     * coordinates as a test of this method.
     ****************************************************/
    std::cout << "\nConverting Galactic to celestial coordinates:" << std::endl;
    std::cout << "   GLon: " << glon << " (degrees)" << std::endl;
    std::cout << "   GLat: " << glat << " (degrees)" << std::endl;
    WcsHandler::Galactic2Celestial(&glon, &glat) ;
    // Now print out the resulting computation
    std::cout << "   Converted coordinates:" << std::endl;
    std::cout << "   RA : " << glon << " (degrees)" << std::endl;
    std::cout << "   Dec: " << glat << " (degrees)" << std::endl;
    
    /****************************************************
     * Test the angular separation calculation
     ****************************************************/
    double RA1(RA), Dec1(Dec), RA2(RA), Dec2(Dec-90.0) ;
    
    std::cout << "\nAngular separation:" << std::endl;
    std::cout << "   RA 1: " << RA1 << std::endl;
    std::cout << "   Dec1: " << Dec1 << std::endl;
    std::cout << "   RA 2: " << RA2 << std::endl;
    std::cout << "   Dec2: " << Dec2 << std::endl;
    std::cout << "   Angular Separation: " << WcsHandler::AngularSeparation_Deg(RA1, Dec1, RA2, Dec2) << std::endl;
    
    return 0 ;
}
