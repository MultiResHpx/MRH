#include "MultiResHpx.h"

// HEALPix includes
#include "healpix_map.h"
#include "pointing.h"
#include "lsconstants.h"
#include "rangeset.h"

// STD includes
#include <time.h>
#include <vector>

void DiscQueryTest(int order,float thetaDeg,float phiDeg,float radius);
void DiscIncQueryTest(int order,float thetaDeg,float phiDeg,float radius);
void TriangleQueryTest(int order,float thetaDeg,float phiDeg);
void TriangleIncQueryTest(int order,float thetaDeg,float phiDeg);
void StripQueryTest(int order,float thetaDeg,float phiDeg);
void StripIncQueryTest(int order,float thetaDeg,float phiDeg);
void NeighborsQueryTest(int order,float thetaDeg,float phiDeg);