#include "MultiResHpx_Map.h"
//#include "MultiResHpx.h"

#include "Common.h"

// HEALPix includes
#include "healpix_map.h"
#include "pointing.h"
#include "lsconstants.h"

// STD includes
#include <time.h>
#include <vector>

void PointQueryTest(int order,float latitude,float longitude);
void DataPointQueryTest(int hpxOrder);
void InterpolationQueryTest(int order,float thetaDeg,float phiDeg);