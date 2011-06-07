cernlib = -L/usr/lib -lkernlib -lmathlib -lpacklib -lgraflib -lgrafX11 -lpawlib

main:: src/main.f src/beamGen.f src/beamDist.f src/dqdt.f src/fields.f src/intd2r.f
	gfortran -Wall src/main.f src/beamGen.f src/beamDist.f src/dqdt.f src/fields.f src/intd2r.f -o fel ${cernlib}

test:: src/fel.f src/beamdist.f
	gfortran -Wall src/fel.f src/beamdist.f -o fel ${cernlib}

intd2r:: src/intd2r.f
	gfortran -Wall src/intd2r.f -o intd2r ${cernlib}

dqdt:: src/dqdt.f
	gfortran -Wall src/dqdt.f -o dqdt ${cernlib}

beamdist:: src/beamgen.f src/beamdist.f
	gfortran -Wall src/beamgen.f src/beamdist.f -o beamdist ${cernlib}

field:: src/fields.f
	gfortran -Wall src/fields.f -o fields ${cernlib}


