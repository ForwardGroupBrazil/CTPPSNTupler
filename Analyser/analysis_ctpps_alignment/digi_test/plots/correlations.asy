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

	xTicksDef = LeftTicks(1., 0.5);
	yTicksDef = RightTicks(1., 0.5);

	TH2_x_min = 0; TH2_x_max = 10;
	TH2_y_min = -3; TH2_y_max = +3;

	real x_cut = 0;
	if (rps[rpi] == "3")
		x_cut = 1.5;
	if (rps[rpi] == "2")
		x_cut = 3.6;
	
	NewPad("$x\ung{mm}$", "$y\ung{mm}$", axesAbove=true);
	scale(Linear, Linear, Log);
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h2_tracks_x_y"), "def");
	draw((x_cut, -3)--(x_cut, +3), magenta+2pt);
	limits((0, -3), (10, +3), Crop);

	NewPad("$x\ung{mm}$", "$y\ung{mm}$", axesAbove=true);
	scale(Linear, Linear, Log);
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h2_digi_corr_x_y"), "def");
	draw((x_cut, -3)--(x_cut, +3), magenta+2pt);
	limits((0, -3), (10, +3), Crop);

	NewPad("$x\ung{mm}$", "$y\ung{mm}$", axesAbove=true);
	scale(Linear, Linear, Log);
	draw(RootGetObject(f, "RP "+rps[rpi]+"/h2_digi_w_corr_x_y"), "def");
	draw((x_cut, -3)--(x_cut, +3), magenta+2pt);
	limits((0, -3), (10, +3), Crop);

	//yaxis(XEquals(exp_V, false), black+2pt+dashed);
}
