import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
//	"period1_alignment/10077",
//	"period1_alignment/10079",
//	"period1_alignment/10080",
//	"period1_alignment/10081",
//	"period1_alignment/10082",

	"period1_physics/fill_5261",
	"period1_physics/fill_5264",
	"period1_physics/fill_5265",
	"period1_physics/fill_5266",
	"period1_physics/fill_5267",
	"period1_physics/fill_5274",
	"period1_physics/fill_5275",
	"period1_physics/fill_5276",
	"period1_physics/fill_5277",
	"period1_physics/fill_5279",
};

string rps[] = {
	"L_1_F",
	"L_1_N",
	"R_1_N",
	"R_1_F",
};

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int rpi : rps.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rps[rpi], "_", "\_") + "}");
}

for (int dsi : datasets.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX " + replace(datasets[dsi], "_", "\_") + "}");

	for (int rpi : rps.keys)
	{
		NewPad("$x\ung{mm}$", "$y\ung{mm}$");

		scale(Linear, Linear, Log);

		TH2_x_min = 5; TH2_x_max = 13;
		TH2_y_min = -3; TH2_y_max = +3;

		TH2_z_min = 1; TH2_z_max = 10^4;

		RootGetObject(topDir + datasets[dsi]+"/distributions.root", "after selection/h2_y_vs_x_" + rps[rpi] + "_sel");
		draw(robj, "def");

		//limits((0, -15), (20, +15), Crop);
	}
}

GShipout(hSkip=1mm, vSkip=1mm);
