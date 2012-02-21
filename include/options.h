#ifndef OPTIONS_H
#define OPTIONS_H

#include <iostream>
#include <cstdlib>

struct Options
{
    std::string inputFileName;
    std::string outputFileName;
    double beamEnergy;
    double bex;
    double bey;
    double emx;
    double emy;
    double lb;
    double lu;
    double bMax;
    double bunchCharge;
    int n;
    int timeStep;
    int mMax;
    int nMax;
    int xGrid;
    int yGrid;
    Options(int, char**, int);
};

Options::Options(int argc, char** argv, int lastargs):
    inputFileName("./input/run2.inp"),
    outputFileName("./data/bunch.txt"),
    beamEnergy(17.5e9),
    bex(32.0),
    bey(32.0),
    emx(1.4e-6),
    emy(1.4e-6),
    lb(25.0),
    lu(35.6),
    bMax(1.0),
    bunchCharge(1.0e-9),
    n(1000),
    timeStep(10),
    mMax(50),
    nMax(50),
    xGrid(1000),
    yGrid(1000)
{
    int limit = argc - lastargs;
    std::string s;
    int i = 1;
    while( i < limit )
    {
        s.assign( argv[i++] );
        if( '-' == s[0] )
        {
            s.assign( s.substr(1) );
            if( s == "ifile" )
            {
                if( i < limit)
                {
                    inputFileName = std::string(argv[i++]);
                }
            }
            else if ( s == "ofile" )
            {
                if( i < limit)
                {
                    outputFileName = std::string(argv[i++]);
                }
            }
            else if ( s == "energy" )
            {
                if( i < limit)
                {
                    beamEnergy = atof(argv[i++]);
                }
            }
            else if ( s == "charge" )
            {
                if( i < limit)
                {
                    bunchCharge = atof(argv[i++]);
                }
            }
            else if ( s == "nparticle" )
            {
                if( i < limit)
                {
                    n = atof(argv[i++]);
                }
            }
            else if ( s == "step" )
            {
                if( i < limit)
                {
                    timeStep = atof(argv[i++]);
                }
            }
            else if ( s == "mmode" )
            {
                if( i < limit)
                {
                    mMax = atof(argv[i++]);
                }
            }
            else if ( s == "nmode" )
            {
                if( i < limit)
                {
                    nMax = atof(argv[i++]);
                }
            }
            else
            {
                std::cerr << "\n*** ERROR *** Unrecognized option: " << s << std::endl;
            }
        }
        else
        {
            std::cerr << "\n*** ERROR *** Unable to interpret command line argument: " << s << std::endl;
        }
    }
}

#endif // OPTIONS_H
