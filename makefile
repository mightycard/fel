cernlib = -L/usr/lib -lkernlib -lmathlib -lpacklib -lgraflib -lgrafX11 -lpawlib
#cern77 = -L/cern/2003/lib -lkernlib -lmathlib -lpacklib -lgraflib -lgrafX11 -lpawlib

test:: src/test.f
	gfortran -Wall src/test.f -o test ${cernlib}

fel:: src/fel.f src/beamdist.f
	g77 -Wall src/fel.f src/beamdist.f -o fel ${cern77}

main:: src/main.f src/beamGen.f src/beamDist.f src/dqdt.f src/fields.f src/intd2r.f
	gfortran -Wall src/main.f src/beamGen.f src/beamDist.f src/dqdt.f src/fields.f src/intd2r.f -o fel ${cernlib}

main2:: src/main.f src/beamgen.f src/beamdist.f src/dqdt.f src/fields.f src/intd2r.f
	gfortran -Wall src/main.f src/beamgen.f src/beamdist.f src/dqdt.f src/fields.f src/intd2r.f -o fel ${cernlib}

intd2r:: src/intd2r.f
	gfortran -Wall src/intd2r.f -o intd2r ${cernlib}

dqdt:: src/dqdt.f
	gfortran -Wall src/dqdt.f -o dqdt ${cernlib}

beamdist:: src/beamgen.f src/beamdist.f
	gfortran -Wall src/beamgen.f src/beamdist.f -o beamdist ${cernlib}

field:: src/sc.f
	gfortran -Wall src/fields.f -o fields ${cernlib}


