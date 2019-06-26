# SLDecayFinder
A Marlin processor to identify semi-leptonic decay of heavy (B/C hadron) flavours in jets and count nSLDecays

# SLDecayFinder is a Marlin processor that looks for semi-leptonic decays (currently at generator level) and gives out the number of semi-leptonic decays of heavy flavour(B/C) hadrons in jets.

# The processor can be followed by additional processors to do analysis on semi-leptonic decays.

# the most frequent use of this processor is to correct neutrino energies.

# The output of the processor is a collection that contains output parameters (number of semi-leptonic decays and index of decaying hadron in MCParticle collection) and is input of the second processor (NuCorrector, for instance)

Quick steps to build the mymarlin example:
------------------------------------------

    . /path/to/ilcsoft/installation/v01-XX/init_ilcsoft.sh
    mkdir build
    cd build
    cmake -C $ILCSOFT/ILCSoft.cmake ..
    make install # this should create library: lib/libmymarlin.so



Quick steps to load mymarlin example into Marlin:
-------------------------------------------------

    export MARLIN_DLL=$PWD/lib/libmymarlin.so
    Marlin -x > mysteer.xml
    grep mymarlin mysteer.xml # should display ... Loading shared library ... (lib/libmymarlin.so)



Quick steps to change this example into your own Marlin plugin:
---------------------------------------------------------------

    1.) rename source files:

        mv include/MyProcessor.h include/ChooseAReasonableNameForThisClass.h
        mv src/MyProcessor.cc src/ChooseAReasonableNameForThisClass.cc

    2.) change CMakeLists.txt:

        PROJECT( ChooseAReasonableName )

        Check DEPENDENCIES for additional required / optional packages


    3.) compile your new plugin and load it into Marlin as described above




Infos and support:
------------------

iLCSoft general documentation:
http://ilcsoft.desy.de/portal/general_documentation/index_eng.html

CMake website:
http://www.cmake.org

The Linear Collider Forum:
http://forum.linearcollider.org
