import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
//	"period1_alignment/10077",

	"period1_physics_margin/fill_4947",
	//"period1_physics_margin/fill_4953",
	"period1_physics/fill_4985",
	"period1_physics/fill_5021",
	"period1_physics/fill_5045",

	"period1_physics/fill_5261",
	"period1_physics/fill_5288",
};

/*
string datasets[] = {
	"period2_physics/fill_5393",
	"period2_physics/fill_5423",
	"period2_physics/fill_5441",
	"period2_physics/fill_5451",
};
*/

string rps[], rp_labels[];
real rp_x_min[], rp_x_max[];

rps.push("3"); rp_labels.push("L-210-fr-hr"); rp_x_min.push(8); rp_x_max.push(18);
rps.push("2"); rp_labels.push("L-210-nr-hr"); rp_x_min.push(7); rp_x_max.push(15);
//rps.push("102"); rp_labels.push("R-210-nr-hr"); rp_x_min.push(7); rp_x_max.push(16.5);
//rps.push("103"); rp_labels.push("R-210-fr-hr"); rp_x_min.push(6.5); rp_x_max.push(16.5);

//yTicksDef = RightTicks(0.5, 0.1);

//----------------------------------------------------------------------------------------------------

//NewPad(false);
for (int rpi : rps.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rp_labels[rpi], "_", "\_") + "}");
}

NewRow();

frame f_leg;

for (int rpi : rps.keys)
{
	NewPad("$x\ung{mm}$");

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		transform t = yscale(1);
		
		/*
		if (qi == 0 && dsi == 0)
			t = yscale(0.5);

		if (qi == 0 && dsi == 2)
			t = yscale(0.38) * shift(0.5, 0);
		*/

		//TH1_x_min = rp_x_min[rpi];
		TH1_x_min = 1;
		//TH1_x_max = rp_x_max[rpi];

		RootGetObject(topDir + datasets[dsi]+"/reconstruction_test.root", "none/without cuts/h_x_" + rps[rpi]);
		draw(t, robj, "vl", p, replace(datasets[dsi], "_", "\_"));
	}

	xlimits(5, 20, Crop);
	//limits((0, 0), (20, +3), Crop);

	f_leg = BuildLegend();
}

//NewRow();
NewPad(false);
add(f_leg);

GShipout(hSkip=1mm, vSkip=1mm);
