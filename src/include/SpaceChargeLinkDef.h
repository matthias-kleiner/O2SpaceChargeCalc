// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SpaceChargeLinkDef_O2.h
/// \author Matthias Kleiner

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedef;

#pragma link C++ class RegularGrid3D < float, 17, 17, 90> + ;
#pragma link C++ class DataContainer3D < float, 17, 17, 90> + ;
#pragma link C++ class O2TPCPoissonSolver < float, 17, 17, 90> + ;

#pragma link C++ class RegularGrid3D < double, 17, 17, 90> + ;
#pragma link C++ class DataContainer3D < double, 17, 17, 90> + ;
#pragma link C++ class O2TPCPoissonSolver < double, 17, 17, 90> + ;

#pragma link C++ class RegularGrid3D < float, 129, 129, 180> + ;
#pragma link C++ class DataContainer3D < float, 129, 129, 180> + ;
#pragma link C++ class O2TPCPoissonSolver < float, 129, 129, 180> + ;

#pragma link C++ class RegularGrid3D < double, 129, 129, 180> + ;
#pragma link C++ class DataContainer3D < double, 129, 129, 180> + ;
#pragma link C++ class O2TPCPoissonSolver < double, 129, 129, 180> + ;

#pragma link C++ class RegularGrid3D < float, 65, 65, 90> + ;
#pragma link C++ class DataContainer3D < float, 65, 65, 90> + ;
#pragma link C++ class O2TPCPoissonSolver < float, 65, 65, 90> + ;

#pragma link C++ class RegularGrid3D < float, 257, 257, 360> + ;
#pragma link C++ class DataContainer3D < float, 257, 257, 360> + ;
#pragma link C++ class O2TPCPoissonSolver < float, 257, 257, 360> + ;

#pragma link C++ class RegularGrid3D < float, 129, 129, 360> + ;
#pragma link C++ class DataContainer3D < float, 129, 129, 360> + ;
#pragma link C++ class O2TPCPoissonSolver < float, 129, 129, 360> + ;

#pragma link C++ class Vector < float, 3> + ;

#endif
