#include "MortonLQT64.h"
// STD includes
#include <time.h>
#include <vector>
#include <math.h>
#include "lsconstants.h"

#define MORTONLQTDEBUG true

void TestHpxToMorton(int64 hpxID, int order);
void TestMortonToHpx(Morton morton);
void TestAllNSIDE8HpxToMorton(int order);
void TestPhiThetaToMorton(double longitude, double latitude, int level);
void TestMortonToPhiTheta(Morton morton);
void TestParentOfMorton(Morton morton);
void TestChildrenOfMorton(Morton morton);
void TestSiblingsOfMorton(Morton morton);
void TestInsertNode(int maxdepth,int num_insert, char* argv[],int isRadians);
void TestSearchNode(int maxdepth,int num_insert, char* argv[]);
void TestDeleteNode(int maxdepth,int num_insert, char* argv[]);
void TestWriteTree(int maxdepth,int num_insert, std::string filename,char* argv[]);
void TestLoadTree(int maxdepth,std::string filename);


void TestMask();
