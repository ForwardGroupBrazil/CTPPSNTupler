#include <map>

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

struct TrackData
{
	bool valid = false;

	// track position, in mm
	double x = 0.;
	double y = 0.;

	// track angle, in rad
	double th_x = 0.;
	double th_y = 0.;

	void operator= (const TotemRPLocalTrack &ft)
	{
		valid = ft.isValid();

		x = ft.getX0();
		y = ft.getY0();

		th_x = ft.getTx();
		th_y = ft.getTy();
	}
};

//----------------------------------------------------------------------------------------------------

struct TrackDataCollection : public map<unsigned int, TrackData>
{
};
