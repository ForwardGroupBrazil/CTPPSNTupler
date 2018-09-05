import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"period1_physics_margin/fill_4947",
	"period1_physics/fill_5261",
};

string rps[], rp_labels[];
rps.push("3"); rp_labels.push("L-210-fr-hr");
rps.push("2"); rp_labels.push("L-210-nr-hr");
rps.push("102"); rp_labels.push("R-210-nr-hr");
rps.push("103"); rp_labels.push("R-210-fr-hr");

string alignments[] = {
	"none",
//	"method x",
};

string selection_dir = "without cuts";

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int rpi : rps.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rp_labels[rpi], "_", "\_") + "}");
}

for (int dsi : datasets.keys)
{
	for (int ai : alignments.keys)
	{
		NewRow();

		NewPad(false);
		label("\vbox{\SetFontSizesXX\hbox{" + replace(datasets[dsi], "_", "\_") + "}\hbox{alignment: " + alignments[ai] + "}}");

		for (int rpi : rps.keys)
		{
			NewPad("$x\ung{mm}$", "$y\ung{mm}$");
	
			scale(Linear, Linear, Log);	
				
			TH2_x_min = 5; TH2_x_max = 10;
			TH2_y_min = -2; TH2_y_max = +2;
	
			TH2_z_min = 1; TH2_z_max = 10^5;
	
			RootGetObject(topDir + datasets[dsi]+"/reconstruction_test.root", alignments[ai] + "/" + selection_dir + "/h2_y_vs_x_" + rps[rpi]);
			draw(robj, "def");
	
			//limits((0, -15), (20, +15), Crop);
		}
	}
}

GShipout(hSkip=1mm, vSkip=1mm);
