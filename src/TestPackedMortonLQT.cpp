#include "TestMortonLQT.h"


#ifdef NOT_IMPLEMENTED

void TestMask()
{
   int max_order = 64;
	uint64 mask = 0;
	cout << sizeof(int64) << " " << sizeof(int) << endl;
	uint64 to_shift = 0;
	for(int pos = 1; pos <= max_order; pos++)
	{
	   //shift_amount = ((2*i)-2);
      //   mask = 3 << shift_amount;
	   //cout << "Mask @ Order: " << i << " Shift Amount: " << shift_amount << " : " << mask << endl;
		if (pos == 1){ to_shift = 1 << pos; }
		 else if (pos == 2){ to_shift = 1 << pos;  }
		 else if (pos == 3){ to_shift = 1 << pos;   }
		 else if (pos == 4){ to_shift = 1 << pos;   }
		 else if (pos == 5){ to_shift = 1 << pos;   }
		 else if (pos == 6){ to_shift = 1 << pos;   }
		 else if (pos == 7){ to_shift = 1 << pos;   }
		 else if (pos == 8){ to_shift = 1 << pos;   }
		 else if (pos == 9){ to_shift = 1 << pos;   }
		 else if (pos == 10){ to_shift = 1 << pos;   }
		 else if (pos == 11){ to_shift = 1 << pos;   }
		 else if (pos == 12){ to_shift = 1 << pos;   }
		 else if (pos == 13){ to_shift = 1 << pos;   }
		 else if (pos == 14){ to_shift = 1 << pos;   }
		 else if (pos == 15){ to_shift = 1 << pos;   }
		 else if (pos == 16){ to_shift = 1 << pos;   }
		 else if (pos == 17){ to_shift = 1 << pos;   }
		 else if (pos == 18){ to_shift = 1 << pos;   }
		 else if (pos == 19){ to_shift = 1 << pos;   }
		 else if (pos == 20){ to_shift = 1 << pos;    }
		 else if (pos == 21){ to_shift = 1 << pos;    }
		 else if (pos == 22){ to_shift = 1 << pos;    }
		 else if (pos == 23){ to_shift = 1 << pos;    }
		 else if (pos == 24){ to_shift = 1 << pos;    }
		 else if (pos == 25){ to_shift = 1 << pos;    }
		 else if (pos == 26){ to_shift = 1 << pos;   }
		 else if (pos == 27){ to_shift = 1 << pos;   }
		 else if (pos == 28){ to_shift = 1 << pos;   }
		 else if (pos == 29){ to_shift = 1 << pos;  }
		 else if (pos == 30){ to_shift = 1 << pos;   }
		 else if (pos == 31){ to_shift = 1 << 1; to_shift = to_shift << 30;   }
		 else if (pos == 32){ to_shift = 1 << 1; to_shift = to_shift << 31;   }
		 else if (pos == 33){ to_shift = 1 << 1; to_shift = to_shift << 32;   }
		 else if (pos == 34){ to_shift = 1 << 2; to_shift = to_shift << 32;   }
		 else if (pos == 35){ to_shift = 1 << 3; to_shift = to_shift << 32;   }
		 else if (pos == 36){ to_shift = 1 << 4; to_shift = to_shift << 32;   }
		 else if (pos == 37){ to_shift = 1 << 5; to_shift = to_shift << 32;   }
		 else if (pos == 38){ to_shift = 1 << 6; to_shift = to_shift << 32;   }
		 else if (pos == 39){ to_shift = 1 << 7; to_shift = to_shift << 32;   }
		 else if (pos == 40){ to_shift = 1 << 8; to_shift = to_shift << 32;   }
		 else if (pos == 41){ to_shift = 1 << 9; to_shift = to_shift << 32;   }
		 else if (pos == 42){ to_shift = 1 << 10; to_shift = to_shift << 32;   }
		 else if (pos == 43){ to_shift = 1 << 11; to_shift = to_shift << 32;   }
		 else if (pos == 44){ to_shift = 1 << 12; to_shift = to_shift << 32;   }
		 else if (pos == 45){ to_shift = 1 << 13; to_shift = to_shift << 32;   }
		 else if (pos == 46){ to_shift = 1 << 14; to_shift = to_shift << 32;   }
		 else if (pos == 47){ to_shift = 1 << 15; to_shift = to_shift << 32;   }
		 else if (pos == 48){ to_shift = 1 << 16; to_shift = to_shift << 32;   }
		 else if (pos == 49){ to_shift = 1 << 17; to_shift = to_shift << 32;   }
		 else if (pos == 50){ to_shift = 1 << 18; to_shift = to_shift << 32;   }
		 else if (pos == 51){ to_shift = 1 << 19; to_shift = to_shift << 32;   }
		 else if (pos == 52){ to_shift = 1 << 20; to_shift = to_shift << 32;   }
		 else if (pos == 53){ to_shift = 1 << 21; to_shift = to_shift << 32;   }
		 else if (pos == 54){ to_shift = 1 << 22; to_shift = to_shift << 32;   }
		 else if (pos == 55){ to_shift = 1 << 23; to_shift = to_shift << 32;   }
	  	 else if (pos == 56){ to_shift = 1 << 24; to_shift = to_shift << 32;   }
		 else if (pos == 57){ to_shift = 1 << 25; to_shift = to_shift << 32;   }
		 else if (pos == 58){ to_shift = 1 << 26; to_shift = to_shift << 32;   }
		 else if (pos == 59){ to_shift = 1 << 27; to_shift = to_shift << 32;   }
		 else if (pos == 60){ to_shift = 1 << 28; to_shift = to_shift << 32;   }
		 else if (pos == 61){ to_shift = 1 << 29; to_shift = to_shift << 32;   }
		 else if (pos == 62){ to_shift = to_shift << 30; to_shift = to_shift << 32;   }
		 else if (pos == 63){ to_shift = to_shift << 31; to_shift = to_shift << 32;   }
		 else if (pos == 64){ to_shift = to_shift << 32; to_shift = to_shift << 32;   }
       //mask = 3 << to_shift;
	    cout << "Mask @ Order: " << pos << " Mask: " << to_shift << endl;
	
	}
}



//
// TEST HpxToMorton, MortonToHpx
//
void TestHpxToPackedMorton(int hpxID, int order)
{
  //Inputs
  // arg[1] = HPX index
  // arg[2] = HPX order 
  // 2-bit address to Morton table
  //
  // MSB  LSB  Morton
  //============
  // 0    0   1
  // 0    1   2
  // 1    0   3
  // 1    1   4
  MortonLQT myMLQT = MortonLQT(13);   
  cout << "#### TEST #1 FROM COMMAND LINE ####" << endl;
  cout << "Inputs:" << endl;
  cout << "  HPX id: " << hpxID << endl;
  cout << "  Order: " << order << endl;
  cout << "Number 2-bit groups: " << order << endl;
  cout << "*** Convert HPX -> Morton ***" << endl;
  PackedMorton morton = myMLQT.HpxToPackedMorton(hpxID,order);                   
  cout << "Morton Code: "; PrintPackedMorton(morton); cout << endl;
}

void TestPackedMortonToHpx(PackedMorton morton)
{
  MortonLQT myMLQT = MortonLQT(13);   
  int64 hpxId;
  int order;
  cout << "#### TEST #2 FROM COMMAND LINE ####" << endl;
  cout << "Inputs:" << endl;
  cout << "  Morton: "; PrintPackedMorton(morton); cout << endl;
  cout << "*** Convert Morton -> HPX ***" << endl;
  myMLQT.PackedMortonToHpx(morton,hpxId,order);
       
  cout << "HPX id: " << hpxId << " order: " << order << endl;
}


void TestLongLatToPackedMorton(float longitude, float latitude, int level)
{     
    MortonLQT myMLQT = MortonLQT(29);   
	cout << "Inputs: " << endl;
	cout << "   Longitude: " << longitude << endl;
	cout << "   Latitude:  " << latitude << endl;
	cout << "   Morton Level: " << level << endl;

	PackedMorton morton = myMLQT.LongLatToPackedMorton(longitude,latitude,level);
	cout << "==================" << endl;
	cout << "   Unpacked Morton Code: "; PrintPackedMorton(morton); cout << endl;
	cout << "SizeOfMorton: " << sizeof(morton);
}

void TestPackedMortonToLongLat(PackedMorton morton)
{    
    MortonLQT myMLQT = MortonLQT(13);   
	float longHpx,coLatHpx;
	cout << "Inputs:" << endl;
	cout << "   Morton: "; PrintPackedMorton(morton); cout << endl;

	myMLQT.PackedMortonToLongLat(morton,longHpx,coLatHpx);
	cout << "==================" << endl;
	cout << "   Longitude: " << longHpx << endl;
	cout << "   Latitude: " << coLatHpx << endl;
}

#endif