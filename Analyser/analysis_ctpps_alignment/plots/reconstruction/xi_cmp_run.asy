import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"period1_alignment/10077",	// reference

	"period1_physics_margin/fill_4947",
//	"period1_physics_margin/fill_4953",
//	"period1_physics_margin/fill_4961",
//	"period1_physics_margin/fill_4964",
//	"period1_physics/fill_4964",
//	"period1_physics_margin/fill_4976",
	"period1_physics/fill_4985",
//	"period1_physics/fill_4988",
//	"period1_physics/fill_4990",

//	"period1_physics/fill_5005",
//	"period1_physics/fill_5013",
	"period1_physics/fill_5017",
//	"period1_physics/fill_5020",
//	"period1_physics/fill_5021",
//	"period1_physics/fill_5024",
//	"period1_physics/fill_5026",
//	"period1_physics/fill_5027",
//	"period1_physics/fill_5028",
//	"period1_physics/fill_5029",
	"period1_physics/fill_5030",
//	"period1_physics/fill_5038",
//	"period1_physics/fill_5043",
//	"period1_physics/fill_5045",
//	"period1_physics/fill_5048",
//	"period1_physics/fill_5052",

//	"period1_physics/fill_5261",
//	"period1_physics/fill_5264",
//	"period1_physics/fill_5265",
//	"period1_physics/fill_5266",	//
//	"period1_physics/fill_5267",
//	"period1_physics/fill_5274",
//	"period1_physics/fill_5275",
//	"period1_physics/fill_5276",
//	"period1_physics/fill_5277",	//
//	"period1_physics/fill_5279",
//	"period1_physics/fill_5287",
//	"period1_physics/fill_5288",

//	"period1_alignment/10077",
//	"period1_alignment/10081",
};

/*
string datasets[] = {
	"period1_alignment/10077",	// reference

	"period1_physics/fill_5043",
	//"period1_physics/fill_5045",
	"period1_physics/fill_5048",
	"period1_physics/fill_5265",
	//"period1_physics/fill_5266",
	"period1_physics/fill_5267",
};
*/

int rp_ids[];
string rp_labels[];
real rp_norm_min[];
real rp_norm_max[];

rp_ids.push(3); rp_labels.push("L-210-fr-hr"); rp_norm_min.push(0.084); rp_norm_max.push(0.114);
rp_ids.push(2); rp_labels.push("L-210-nr-hr"); rp_norm_min.push(0.083); rp_norm_max.push(0.113);
rp_ids.push(102); rp_labels.push("R-210-nr-hr"); rp_norm_min.push(0.085); rp_norm_max.push(0.140);
rp_ids.push(103); rp_labels.push("R-210-fr-hr"); rp_norm_min.push(0.085); rp_norm_max.push(0.145);

string alignments[] = {
//	"none",
	"method x",
//	"method y",
};

string cut_option = "without cuts";

bool cropToDetails = false;

real x_min[] = {0.05, 0.05, 0.05,  0.06};
real x_max[] = {0.11, 0.11, 0.14,  0.14};
real y_min[] = {0.03, 0.03, 0.015, 0.015};
real y_max[] = {0.04, 0.04, 0.020, 0.020};


xSizeDef = 8cm;
xTicksDef = (cropToDetails) ? LeftTicks(0.02, 0.01) : LeftTicks(0.05, 0.01);

//----------------------------------------------------------------------------------------------------

real GetNormalisation(RootObject obj, real xi_min, real xi_max)
{
	RootObject x_axis = obj.oExec("GetXaxis");
	int bin_min = x_axis.iExec("FindBin", xi_min);
	int bin_max = x_axis.iExec("FindBin", xi_max);

	real S = 0.;
	for (int bin = bin_min; bin <= bin_max; ++bin)
	{
		S += obj.rExec("GetBinContent", bin);
	}

	return S;
}

//----------------------------------------------------------------------------------------------------

/*
NewPad(false);
label(cut_option);

NewRow();

NewPad(false);
for (int rpi : rp_ids.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rp_labels[rpi], "_", "\_") + "}");
}
*/

for (int ai : alignments.keys)
{
	NewRow();

	//NewPad(false);
	//label("{\SetFontSizesXX " + alignments[ai] + "}");

	frame f_legend;

	for (int rpi : rp_ids.keys)
	{
		if (rpi == 2)
			NewRow();

		NewPad("$\xi$");

		for (int dsi : datasets.keys)
		{
			bool alignmentRun = (find(datasets[dsi], "period1_alignment") != -1);

			string alignments_eff = alignments[ai];
			if (alignmentRun)
				alignments_eff = "none";

			string f = topDir + datasets[dsi] + "/reconstruction_test.root";
			RootObject obj = RootGetObject(f, alignments_eff+"/"+cut_option+"/" + format("h_xi_%u", rp_ids[rpi]));

			real norm = GetNormalisation(obj, rp_norm_min[rpi], rp_norm_max[rpi]);

			pen p = (alignmentRun) ? black+2pt : StdPen(dsi);

			draw(scale(1., 1./norm), obj, "vl", p, replace(datasets[dsi], "_", "\_"));
		}

		if (cropToDetails)
			limits((x_min[rpi], y_min[rpi]), (x_max[rpi], y_max[rpi]), Crop);

		f_legend = BuildLegend();

		currentpicture.legend.delete();
		AttachLegend(rp_labels[rpi]);
	}

	NewPad(false);
	attach(f_legend);
}


GShipout(hSkip=1mm, vSkip=1mm, margin=1mm);
