import root;
import pad_layout;

string rps[];
string rp_labels[];

rps.push("3"); rp_labels.push("L, 210, fr, hr");
rps.push("2"); rp_labels.push("L, 210, nr, hr");

string plane_U = "1";
string plane_V = "0";

string f = "../distributions.root";

xSizeDef = 8cm;

TH2_palette = Gradient(white, blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

for (int rpi : rps.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX RP " + rp_labels[rpi] + "}");

	real exp_X = 0, exp_Y = 0;
	real rod_U_x = 0, rod_U_y;
	real rod_V_x = 0, rod_V_y;
	if (rps[rpi] == "3")
	{
		exp_X = 3.35; exp_Y = -0.1;
		rod_U_x = 0.713; rod_U_y = 0.702;
		rod_V_x = 0.702; rod_U_y = -0.713;
	}
	if (rps[rpi] == "2")
	{
		exp_X = 4.0; exp_Y = -0.1;
		rod_U_x = 0.624; rod_U_y = 0.781;
		rod_V_x = 0.781; rod_U_y = -0.624;
	}

	real exp_U = rod_U_x*exp_X + rod_U_x*exp_Y;
	real exp_V = rod_V_x*exp_X + rod_V_x*exp_Y;

	xTicksDef = LeftTicks(1., 0.5);
	yTicksDef = RightTicks(1., 0.5);
	
	NewPad("$x\ung{mm}$", "$y\ung{mm}$", axesAbove=true);
	TH2_x_min = 0; TH2_x_max = 10;
	TH2_y_min = -3; TH2_y_max = +3;
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h2_tracks_x_y"), "def");
	draw((exp_X, exp_Y), mCr+5pt+(black+2pt));
	limits((0, -3), (10, +3), Crop);

	xTicksDef = LeftTicks(5., 1.);
	yTicksDef = RightTicks(200., 100.);

	NewPad("$U\ung{mm}$");
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h_plane_"+plane_U+"_digi_abs"), "vl", blue, "digi (plane " + plane_U + ")");
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h_tracks_U"), "vl", red, "tracks");
	xlimits(-5, 15, Crop);
	AttachLegend("U projection");

	yaxis(XEquals(exp_U, false), black+2pt+dashed);

	NewPad("$V\ung{mm}$");
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h_plane_"+plane_V+"_digi_abs"), "vl", blue, "digi (plane " + plane_V + ")");
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h_tracks_V"), "vl", red, "tracks");
	xlimits(-5, 15, Crop);
	AttachLegend("V projection");

	yaxis(XEquals(exp_V, false), black+2pt+dashed);
}
