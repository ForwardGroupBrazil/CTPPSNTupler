#include <map>
#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct SelectionRange
{
	double min;
	double max;

	SelectionRange(double _mi=0., double _ma=0.) : min(_mi), max(_ma)
	{
	}
};

//----------------------------------------------------------------------------------------------------

// parameter defaults

unsigned int lhc_fill = 0;

bool alignment_run = false;

double n_si = 3.;

bool cut1_apply = true;
bool cut2_apply = true;
bool cut3_apply = false;
bool cut4_apply = false;

double cut1_a = -1., cut1_c = 0., cut1_si = 0.20;
double cut2_a = -1., cut2_c = 0., cut2_si = 0.20;
double cut3_a = -1., cut3_c = 0., cut3_si = 0.20;
double cut4_a = -1., cut4_c = 0., cut4_si = 0.20;

// x range for distribution matching
std::map<std::string, SelectionRange> selectionRangesX;

// x range for y alignment (fit and extrapolation)
std::map<std::string, SelectionRange> yAlignmentRangeX;
