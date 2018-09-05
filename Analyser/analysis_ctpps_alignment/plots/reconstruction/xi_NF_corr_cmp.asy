import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"period1_alignment/10077",
//	"period1_alignment/10081",

	"period1_physics_margin/fill_4947",
//	"period1_physics_margin/fill_4953",
//	"period1_physics_margin/fill_4961",
//	"period1_physics_margin/fill_4964",
	
//	"period1_physics/fill_4964",
//	"period1_physics_margin/fill_4976",

//	"period1_physics/fill_4985",
//	"period1_physics/fill_4988",
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
//	"period1_physics/fill_5048",
	"period1_physics/fill_5052",

	"period1_physics/fill_5261",
//	"period1_physics/fill_5264",
//	"period1_physics/fill_5265",
//	"period1_physics/fill_5266",
//	"period1_physics/fill_5267",

	"period1_physics/fill_5274",
	"period1_physics/fill_5279",
};

//string datasets[] = {
//	"period1_physics_margin/fill_4953",
//};

string datasets[] = {
	"period1_physics_margin/fill_4947",
	"period1_physics/fill_5261",
};

string arm_tags[];
string arm_labels[];

arm_tags.push("L"); arm_labels.push("left arm");
arm_tags.push("R"); arm_labels.push("right arm");


string alignment = "method x";

string cut_option = "with cuts";

xSizeDef = 8cm;
ySizeDef = 8cm;

real xi_min = 0.03;
real xi_max = 0.20;

TH2_palette = Gradient(blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

NewPad(false);
label(cut_option + ", " + alignment);

NewRow();

NewPad(false);
for (int ai : arm_tags.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(arm_labels[ai], "_", "\_") + "}");
}

for (int dsi : datasets.keys)
{
	NewRow();

	NewPad(false);
	string bits[] = split(replace(datasets[dsi], "_", "\_"), "/");
	label("\vbox{\SetFontSizesXX\hbox{" + bits[0] + "}\hbox{" + bits[1] + "}}");

	for (int ai : arm_tags.keys)
	{
		NewPad("$\xi^{\rm N}$", "$\xi^{\rm F}$");

		string arm = arm_tags[ai];
		string f = topDir + datasets[dsi] + "/reconstruction_test.root";
		RootObject obj = RootGetObject(f, alignment+"/"+cut_option+"/h2_xi_" + arm + "_F_vs_xi_" + arm + "_N");

		TH2_x_min = xi_min; TH2_x_max = xi_max;
		TH2_y_min = xi_min; TH2_y_max = xi_max;

		draw(obj);

		draw((xi_min, xi_min)--(xi_max, xi_max), black+2pt+dashed);
	
		limits((xi_min, xi_min), (xi_max, xi_max), Crop);
	}
}

GShipout(hSkip=1mm, vSkip=1mm, margin=1mm);
