import root;
import pad_layout;

string topDir = "../../";

string dataset = "period1_physics/fill_4964";
//string dataset = "period1_physics/fill_5024";
//string dataset = "period1_physics/fill_5052";

string reference = "10081";

string rp = "L_1_F";

xSizeDef = 8cm;
ySizeDef = 6cm;

//----------------------------------------------------------------------------------------------------

/*
NewPad(false);
label(replace(dataset, "_", "\_"));

NewRow();
*/

NewPad(false);
label("{\SetFontSizesXX method x}");

NewPad(false);
label("{\SetFontSizesXX method y}");

NewRow();

NewPad("$x\ung{mm}$", "entries (scaled)");
//currentpad.yTicks = RightTicks(0.5, 0.1);

string p_base = rp + "/" + reference + "/method x/c_cmp|";
draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_ref_sel"), "d0,eb", black);
draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_bef"), "d0,eb", blue);
draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_aft"), "d0,eb", red);

xlimits(4, 17, Crop);


NewPad("$x\ung{mm}$", "std.~dev.~of $y\ung{mm}$");
currentpad.yTicks = RightTicks(0.5, 0.1);

string p_base = rp + "/" + reference + "/method y/c_cmp|";
draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_ref_sel"), "d0,eb", black);
draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_bef"), "d0,eb", blue);
draw(RootGetObject(topDir + dataset+"/match.root", p_base + "h_test_aft"), "d0,eb", red);

limits((4, 0), (17, 3), Crop);


GShipout(hSkip=10mm, vSkip=1mm);
