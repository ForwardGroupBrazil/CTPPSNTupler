#include "common.h"

#include "shared_track.h"
#include "shared_alignment.h"

#include "parameters.h"

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

struct Correction
{
	double corr_102, corr_103;

	void Print() const
	{
		printf("102: %.3f\n", corr_102);
		printf("103: %.3f\n", corr_103);
	}
};

map<unsigned int, Correction> corrections;

void AddAlCorrection(unsigned int fill, double c102, double c103)
{
	corrections[fill] = {c102, c103};
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	// apply settings
	ApplySettings();

	// define corrections
	AddAlCorrection(5261, 0.05, 0.05);
	AddAlCorrection(5265, 0.15, 0.10);
	AddAlCorrection(5266, 0.08, 0.00);
	AddAlCorrection(5267, 0.22, 0.12);

	AddAlCorrection(5274, 0.12, 0.10);
	AddAlCorrection(5275, 0.20, 0.10);
	AddAlCorrection(5276, 0.18, 0.17);
	AddAlCorrection(5277, 0.15, 0.05);

	AddAlCorrection(5279, 0.15, 0.10);
	AddAlCorrection(5287, 0.15, 0.15);
	AddAlCorrection(5288, 0.28, 0.20);

	// read in alignments
	AlignmentResultsCollection alignments;
	alignments.Load("match.out");

	// process alignments
	auto it = corrections.find(lhc_fill);
	if (it == corrections.end())
	{
		printf("No correction to be applied.\n");
	} else {
		printf("Applying correction:\n");
		it->second.Print();

		for (auto &ap : alignments)
		{
			for (auto &rp : ap.second)
			{
				if (rp.first == 102)
					rp.second.sh_x += it->second.corr_102;

				if (rp.first == 103)
					rp.second.sh_x += it->second.corr_103;
			}
		}
	}
	
	// save alignments
	alignments.Write("process_alignments.out");

	return 0;
}
