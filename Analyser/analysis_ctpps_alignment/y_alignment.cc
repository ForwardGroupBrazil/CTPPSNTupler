#include "common.h"

#include "shared_track.h"
#include "shared_alignment.h"

#include "parameters.h"

/*
#include "stat.h"

#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
*/
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

TF1 *ff_pol1 = new TF1("ff_pol1", "[0] + [1]*x");

//----------------------------------------------------------------------------------------------------

void DoYAlignment(TProfile *p, const SelectionRange &range, double &result, double &unc)
{
	printf("range: min = %.2f, max = %.2f\n", range.min, range.max);

	// do main fit
	ff_pol1->SetParameters(0., 0.);
	p->Fit(ff_pol1, "", "", range.min, range.max);

	result = ff_pol1->GetParameter(0);
	double fit_unc = ff_pol1->GetParError(0);

	double slope = ff_pol1->GetParameter(1);

	// save main fit
	p->Write("p_y_vs_x");

	// secondary fit
	ff_pol1->SetParameters(0., 0.);
	p->Fit(ff_pol1, "", "", range.min+1, range.max-1);

	// determine uncertainty
	double x_unc = 150E-3;	// in mm
	double x_scale_unc = slope * x_unc;
	double range_unc = fabs(result - ff_pol1->GetParameter(0));

	unc = sqrt(fit_unc*fit_unc + x_scale_unc*x_scale_unc + range_unc*range_unc);

	// print results
	printf("y(x=0) = %.3f mm\n", result);
	printf("    fit uncertainty = %.3f mm\n", fit_unc);
	printf("    uncertainty from x scale = %.3f mm\n", x_scale_unc);
	printf("    uncertainty from fit range = %.3f mm\n", range_unc);
	printf("    total uncertainty = %.3f mm\n", unc);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	// apply settings
	ApplySettings();

	// list of RPs and their settings
	struct RPData
	{
		string name;
		unsigned int id;
	};

	vector<RPData> rpData = {
		{ "L_1_F", 3 },
		{ "L_1_N", 2 },
		{ "R_1_N", 102 },
		{ "R_1_F", 103 },
	};

	// get input
	TFile *f_in = new TFile("reconstruction_test.root");

	// ouput file
	TFile *f_out = new TFile("y_alignment.root", "recreate");

	// prepare results
	AlignmentResultsCollection results;

	// processing per rp
	for (const auto &rpd : rpData)
	{
		printf("-------------------- %s --------------------\n", rpd.name.c_str());

		TDirectory *rp_dir = f_out->mkdir(rpd.name.c_str());
		gDirectory = rp_dir;

		char buf[100];
		sprintf(buf, "method x/with cuts/p_y_vs_x_%i", rpd.id);
		TProfile *p_input = (TProfile *) f_in->Get(buf);

		auto it_range = yAlignmentRangeX.find(rpd.name);
		if (it_range == yAlignmentRangeX.end())
		{
			printf("ERROR: can't find RP in yAlignmentRangeX.\n");
			continue;
		}

		double result = 0;
		double unc = 0;
		DoYAlignment(p_input, it_range->second, result, unc);

		results["y alignment"][rpd.id] = AlignmentResult(0., 0., result, unc);
	}

	// write results
	FILE *f_results = fopen("y_alignment.out", "w"); 
	results.Write(f_results);
	fclose(f_results);

	// clean up
	delete f_out;
	return 0;
}
