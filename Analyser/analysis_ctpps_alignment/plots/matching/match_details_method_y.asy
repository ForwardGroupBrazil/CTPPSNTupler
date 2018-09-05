import root;
import pad_layout;

string topDir = "../../";

string dataset = "period1_physics/fill_4964";
//string dataset = "period1_physics/fill_5024";
//string dataset = "period1_physics/fill_5052";

string reference = "10081";

string rps[] = {
	"L_1_F",
	"L_1_N",
	"R_1_N",
	"R_1_F",
};

ySizeDef = 5cm;

//----------------------------------------------------------------------------------------------------

NewPad(false);
label(replace(dataset, "_", "\_"));

NewRow();

for (int rpi : rps.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rps[rpi], "_", "\_") + "}");
}

NewRow();

for (int rpi : rps.keys)
{
	NewPad("$x\ung{mm}$", "std.~dev.~of $y\ung{mm}$");
	currentpad.yTicks = RightTicks(0.5, 0.1);

	string p_base = rps[rpi] + "/" + reference + "/method y/c_cmp|";
	draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_ref_sel"), "d0,eb", black);
	draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_bef"), "d0,eb", blue);
	draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_aft"), "d0,eb", red);

	limits((2, 0), (15, 3), Crop);
}

NewRow();

xTicksDef = LeftTicks(1., 0.2);

for (int rpi : rps.keys)
{
	NewPad("shift$\ung{mm}$", "bins in overlap");
	
	string p_base = rps[rpi] + "/" + reference + "/method y";

	RootGetObject(topDir + dataset+"/match.root", p_base + "/g_results");
	real ax[] = { 0. };
	real ay[] = { 0. };
	robj.vExec("GetPoint", 0, ax, ay); real sh_best = ay[0];
	robj.vExec("GetPoint", 1, ax, ay); real sh_best_unc = ay[0];

	draw(RootGetObject(topDir + dataset+"/match.root", p_base + "/g_n_bins"), "p", magenta, mCi+1pt+magenta);

	limits((-5, 0), (+0, 80), Crop);
	yaxis(XEquals(sh_best - sh_best_unc, false), dashed);
	yaxis(XEquals(sh_best, false), solid);
	yaxis(XEquals(sh_best + sh_best_unc, false), dashed);
}

NewRow();

for (int rpi : rps.keys)
{
	NewPad("shift$\ung{mm}$", "$S^2 / N$");
	
	string p_base = rps[rpi] + "/" + reference + "/method y";

	RootGetObject(topDir + dataset+"/match.root", p_base + "/g_results");
	real ax[] = { 0. };
	real ay[] = { 0. };
	robj.vExec("GetPoint", 0, ax, ay); real sh_best = ay[0];
	robj.vExec("GetPoint", 1, ax, ay); real sh_best_unc = ay[0];

	draw(RootGetObject(topDir + dataset+"/match.root", p_base + "/g_chi_sq_norm"), "p", heavygreen, mCi+1pt+heavygreen);

	limits((-5, 0), (+0, 200), Crop);
	yaxis(XEquals(sh_best - sh_best_unc, false), dashed);
	yaxis(XEquals(sh_best, false), solid);
	yaxis(XEquals(sh_best + sh_best_unc, false), dashed);
}

GShipout(hSkip=1mm, vSkip=1mm);
