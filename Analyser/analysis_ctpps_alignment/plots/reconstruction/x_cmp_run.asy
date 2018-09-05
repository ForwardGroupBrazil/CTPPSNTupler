import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"period1_physics_margin/fill_4947",
//	"period1_physics_margin/fill_4953",
//	"period1_physics_margin/fill_4961",
	"period1_physics_margin/fill_4964",
	
//	"period1_physics/fill_4964",
//	"period1_physics_margin/fill_4976",

//	"period1_physics/fill_4985",
	"period1_physics/fill_4988",
//	"period1_physics/fill_4990",
//	"period1_physics/fill_5005",
//	"period1_physics/fill_5013",
//	"period1_physics/fill_5017",
//	"period1_physics/fill_5020",
//	"period1_physics/fill_5021",
//	"period1_physics/fill_5024",
	"period1_physics/fill_5026",
//	"period1_physics/fill_5027",
//	"period1_physics/fill_5028",
//	"period1_physics/fill_5029",
//	"period1_physics/fill_5030",
//	"period1_physics/fill_5038",
//	"period1_physics/fill_5043",
//	"period1_physics/fill_5045",
	"period1_physics/fill_5048",
//	"period1_physics/fill_5052",

	"period1_physics/fill_5261",
//	"period1_physics/fill_5264",
	"period1_physics/fill_5265",
	"period1_physics/fill_5266",
	"period1_physics/fill_5267",

	"period1_alignment/10077",
//	"period1_alignment/10081",
};

string datasets[] = {
	"period1_physics_margin/fill_4947",
//	"period1_physics_margin/fill_4964",

	"period1_physics/fill_4988",
	"period1_physics/fill_5026",
	"period1_physics/fill_5048",

	"period1_physics/fill_5261",
//	"period1_physics/fill_5266",

	"period1_alignment/10077",
//	"period1_alignment/10081",
};

int rp_ids[];
string rp_labels[];
real rp_norm_min[];
real rp_norm_max[];

rp_ids.push(3); rp_labels.push("L-210-fr-hr"); rp_norm_min.push(9.0); rp_norm_max.push(11.0);
rp_ids.push(2); rp_labels.push("L-210-nr-hr"); rp_norm_min.push(9.0); rp_norm_max.push(11.0);
rp_ids.push(102); rp_labels.push("R-210-nr-hr"); rp_norm_min.push(7.0); rp_norm_max.push(10.0);
rp_ids.push(103); rp_labels.push("R-210-fr-hr"); rp_norm_min.push(7.0); rp_norm_max.push(10.0);

string alignments[] = {
//	"none",
	"method x",
//	"method y",
};

string cut_option = "without cuts";

xSizeDef = 8cm;
xTicksDef = LeftTicks(1., 0.5);

//----------------------------------------------------------------------------------------------------

real GetNormalisation(RootObject obj, real x_min, real x_max)
{
	RootObject x_axis = obj.oExec("GetXaxis");
	int bin_min = x_axis.iExec("FindBin", x_min);
	int bin_max = x_axis.iExec("FindBin", x_max);

	real S = 0.;
	for (int bin = bin_min; bin <= bin_max; ++bin)
	{
		S += obj.rExec("GetBinContent", bin);
	}

	return S;
}

//----------------------------------------------------------------------------------------------------

NewPad(false);
label(cut_option);

NewRow();

NewPad(false);
for (int rpi : rp_ids.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rp_labels[rpi], "_", "\_") + "}");
}

for (int ai : alignments.keys)
{
	NewRow();

	NewPad(false);
	label(rotate(90)*Label("{\SetFontSizesXX " + alignments[ai] + "}"));

	for (int rpi : rp_ids.keys)
	{
		NewPad("$x\ung{mm}$");

		TH1_x_min = 1;

		for (int dsi : datasets.keys)
		{
			bool alignmentRun = (find(datasets[dsi], "period1_alignment") != -1);

			string alignments_eff = alignments[ai];
			if (alignmentRun)
				alignments_eff = "none";

			string f = topDir + datasets[dsi] + "/reconstruction_test.root";
			RootObject obj = RootGetObject(f, alignments_eff+"/"+cut_option+"/" + format("h_x_%u", rp_ids[rpi]));

			real norm = GetNormalisation(obj, rp_norm_min[rpi], rp_norm_max[rpi]);

			pen p = (alignmentRun) ? black+2pt : StdPen(dsi+1);

			draw(scale(1., 1./norm), obj, "vl", p, replace(datasets[dsi], "_", "\_"));
		}
	
		//limits((1, 0), (15., 0.10), Crop);
		xlimits(1, 15., Crop);
	}

	frame f_legend = BuildLegend();

	NewPad(false);
	attach(f_legend);
}


GShipout(hSkip=1mm, vSkip=1mm, margin=1mm);
