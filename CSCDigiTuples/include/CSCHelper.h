#ifndef CSCUCLA_CSCDIGITUPLES_CSCHELPER_H
#define CSCUCLA_CSCDIGITUPLES_CSCHELPER_H

//#include "DataFormats/MuonDetId/interface/CSCDetId.h"




namespace CSCHelper {


//Conversion function to check validity
template<typename Target, typename Source>
Target convertTo(Source source, const char name[], bool lenient = false, bool* good = 0)
{
  Target            converted = static_cast<Target>(source);

  if (static_cast<Source>(converted) != source) {
    const Target    lowest    = !std::numeric_limits<Target>::is_signed
                              ? 0
                              : std::numeric_limits<Target>::has_infinity
                              ? -std::numeric_limits<Target>::infinity()
                              :  std::numeric_limits<Target>::min()
                              ;

    std::string problem = "convertTo(): Source value " + std::to_string((double)  source) + " outside of target range "
                         +"["+std::to_string((double)  lowest)+","+std::to_string((double)  std::numeric_limits<Target>::max())+"]"
                         + " for '" + name +"'";

    if (good)      *good      = false;
    if (lenient) {
      std::cerr << "WARNING: " << problem << std::endl;
      return  ( source > static_cast<Source>(std::numeric_limits<Target>::max())
              ? std::numeric_limits<Target>::max()
              : lowest
              );
    }
    throw std::range_error( problem);
  }

  return converted;
}

/*
unsigned short int chamberSerial( CSCDetId id ) {
    int st = id.station();
    int ri = id.ring();
    int ch = id.chamber();
    int ec = id.endcap();
    int kSerial = ch;
    if (st == 1 && ri == 1) kSerial = ch;
    if (st == 1 && ri == 2) kSerial = ch + 36;
    if (st == 1 && ri == 3) kSerial = ch + 72;
    if (st == 1 && ri == 4) kSerial = ch;
    if (st == 2 && ri == 1) kSerial = ch + 108;
    if (st == 2 && ri == 2) kSerial = ch + 126;
    if (st == 3 && ri == 1) kSerial = ch + 162;
    if (st == 3 && ri == 2) kSerial = ch + 180;
    if (st == 4 && ri == 1) kSerial = ch + 216;
    if (st == 4 && ri == 2) kSerial = ch + 234;  // one day...
    if (ec == 2) kSerial = kSerial + 300;
    return convertTo<unsigned short int>(kSerial,"chamberSerial");
}
*/


//written to remove ambiguity between ME11A and ME11B
unsigned int serialize(unsigned int st, unsigned int ri, unsigned int ch, unsigned int ec){
	//set them to be 0 based
	st--;ri--;ch--;ec--;
    //sanity check
    if(ec >= 2 || st >= 4 || ri >= 4 || ch >= 36) {
    	std::cout << "Error: Trying to serialize invalid chamber:" <<
    			" ST = " << st <<
    			" RI = " << ri <<
				" CH = " << ch <<
				" EC = " << ec << std::endl;
    	return 0xffffffff;
    }

	/* EC = 1 -> 2  [1 bit]
	 * ST = 1 -> 4  [2 bits]
	 * RI = 1 -> 4  [2 bits]
	 * CH = 1 -> 36 [6 bits]
	 *
	 * -> need 1+2+2+6 = 11 bits = 2 bytes -> unsigned int
	 *
	 * least significant is channel, goes up to endcap as most significant
	 */
    return  (ec << 10)|(st << 8)|(ri << 6)| ch;
}

struct ChamberId {
	unsigned int station;
	unsigned int ring;
	unsigned int chamber;
	unsigned int endcap;
};

unsigned int serialize(ChamberId id){
    unsigned int st = id.station;
    unsigned int ri = id.ring;
    unsigned int ch = id.chamber;
    unsigned int ec = id.endcap;
    return CSCHelper::serialize(st,ri,ch,ec);
}

/* @brief Unserializes chamber id. Everything <= 1
 *
 */
ChamberId unserialize(unsigned int serial){
	ChamberId c;
	c.chamber		= (serial & 0x0000003f)+1; //6 bits
	serial 			= serial >> 6;
	c.ring 			= (serial & 0x00000003)+1; //2 bits
	serial			= serial >> 2;
	c.station 		= (serial & 0x00000003)+1; //2 bits
	serial			= serial >> 2;
	c.endcap		= serial+1;
	return c;

}


};

#endif /*CSCUCLA_CSCDIGITUPLES_FILLCSCINFO_H*/
