import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
//	"period1_physics_margin/274199",
//	"period1_physics_margin/274241",
//	"period1_physics/274244",
//	"period1_physics/274388",
//	"period1_physics/274958",
//	"period1_physics/274969",
//	"period1_physics/275125",
//	"period1_physics/275310",
//	"period1_physics/275376",
//	"period1_physics/275836",

	"period1_physics_margin/fill_4947",
	"period1_physics_margin/fill_4953",
	"period1_physics_margin/fill_4961",
	"period1_physics_margin/fill_4964",
	
	"period1_physics/fill_4964",
	"period1_physics_margin/fill_4976",

	"period1_physics/fill_4985",
	"period1_physics/fill_4988",
	"period1_physics/fill_4990",
	"period1_physics/fill_5005",
	"period1_physics/fill_5013",
	"period1_physics/fill_5017",
	"period1_physics/fill_5020",
	"period1_physics/fill_5021",
	"period1_physics/fill_5024",
	"period1_physics/fill_5026",
	"period1_physics/fill_5027",
	"period1_physics/fill_5028",
	"period1_physics/fill_5029",
	"period1_physics/fill_5030",
	"period1_physics/fill_5038",
	"period1_physics/fill_5043",
	"period1_physics/fill_5045",
	"period1_physics/fill_5048",
	"period1_physics/fill_5052",
};


string reference = "10081";

string rps[] = {
	"L_1_F",
	"L_1_N",
	"R_1_N",
	"R_1_F",
};

xSizeDef = 8cm;

yTicksDef = RightTicks(0.005, 0.001);

//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int rpi : rps.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rps[rpi], "_", "\_") + "}");
}

for (int dsi : datasets.keys)
{
	string dataset = datasets[dsi];

	NewRow();

	NewPad(false);
	string bits[] = split(replace(dataset, "_", "\_"), "/");
	label("\vbox{\SetFontSizesXX\hbox{" + bits[0] + "}\hbox{" + bits[1] + "}}");
	
	for (int rpi : rps.keys)
	{
		NewPad("$x\ung{mm}$", "entries (scaled)");
	
		string p_base = rps[rpi] + "/" + reference + "/method x/c_cmp|";
		draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_ref_sel"), "d0,eb", black);
		draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_bef"), "d0,eb", blue);
		draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_aft"), "d0,eb", red);
	
		limits((3, 0.005), (17, 0.045), Crop);
	}
}

GShipout(hSkip=1mm, vSkip=1mm);
