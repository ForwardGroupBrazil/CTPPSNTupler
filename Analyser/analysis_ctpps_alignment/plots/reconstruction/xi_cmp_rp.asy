import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"period1_alignment/10077",
	"period1_alignment/10081",

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
};

string datasets[] = {
	"period1_alignment/10077",
//	"period1_physics_margin/fill_4947",
};

int rp_ids[];
string rp_labels[];
pen rp_pens[];

rp_ids.push(3); rp_labels.push("L-210-fr-hr"); rp_pens.push(black);
rp_ids.push(2); rp_labels.push("L-210-nr-hr"); rp_pens.push(red);
rp_ids.push(102); rp_labels.push("R-210-nr-hr"); rp_pens.push(blue);
rp_ids.push(103); rp_labels.push("R-210-fr-hr"); rp_pens.push(heavygreen);

string alignments[] = {
//	"none",
	"method x",
//	"method y",
};

string cut_option = "without cuts";

xSizeDef = 10cm;
xTicksDef = LeftTicks(0.05, 0.01);

//----------------------------------------------------------------------------------------------------

NewPad(false);
label(cut_option);

NewRow();

NewPad(false);
for (int ai : alignments.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(alignments[ai], "_", "\_") + "}");
}

for (int dsi : datasets.keys)
{
	NewRow();

	NewPad(false);
	string bits[] = split(replace(datasets[dsi], "_", "\_"), "/");
	label("\vbox{\SetFontSizesXX\hbox{" + bits[0] + "}\hbox{" + bits[1] + "}}");

	for (int ai : alignments.keys)
	{
		NewPad("$\xi$");
		//scale(Linear, Log(auto));

		for (int rpi : rp_ids.keys)
		{
			string f = topDir + datasets[dsi] + "/reconstruction_test.root";
			RootObject obj = RootGetObject(f, alignments[ai]+"/"+cut_option+"/" + format("h_xi_%u", rp_ids[rpi]));
			draw(obj, "vl", rp_pens[rpi], replace(rp_labels[rpi], "_", "\_"));
		}

		xlimits(0, 0.2, Crop);
		//ylimits(1e1, 2e5, Crop);
	}

	frame f_legend = BuildLegend();

	NewPad(false);
	attach(f_legend);
}


GShipout(hSkip=1mm, vSkip=1mm, margin=1mm);
