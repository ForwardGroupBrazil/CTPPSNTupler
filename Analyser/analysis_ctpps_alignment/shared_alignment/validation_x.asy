import root;
import pad_layout;

string topDir = "../";

string datasets[] = {
	"period1_physics_margin/fill_4947",
	"period1_physics_margin/fill_4953",
	"period1_physics_margin/fill_4961",
	"period1_physics_margin/fill_4964",
	"period1_physics_margin/fill_4976",

	"period1_physics/fill_4964",
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
	"period1_physics/fill_5287",
	"period1_physics/fill_5288",
};

int rp_ids[];
string rps[], rp_labels[];
rp_ids.push(3); rps.push("L_1_F"); rp_labels.push("L-210-fr-hr");
rp_ids.push(2); rps.push("L_1_N"); rp_labels.push("L-210-nr-hr");
rp_ids.push(102); rps.push("R_1_N"); rp_labels.push("R-210-nr-hr");
rp_ids.push(103); rps.push("R_1_F"); rp_labels.push("R-210-fr-hr");

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewRow();

	string dataset = datasets[dsi];

	string f_ref = topDir+dataset+"/reconstruction_test.root";
	string f_test = topDir+dataset+"/shared_reconstruction_test.root";

	NewPad(false);
	label(replace(dataset, "_", "\_"));

	for (int rpi : rps.keys)
	{
		NewPad();

		string rpIdStr = format("%i", rp_ids[rpi]);

		string p_ref = "method x/with cuts/h_x_" + rpIdStr;
		string p_test = "with cuts/h_x_" + rpIdStr;

		RootObject obj_ref = RootGetObject(f_ref, p_ref);
		RootObject obj_test = RootGetObject(f_test, p_test);
		
		draw(obj_ref, "eb", blue);
		draw(obj_test, "eb", red);
	}
}
