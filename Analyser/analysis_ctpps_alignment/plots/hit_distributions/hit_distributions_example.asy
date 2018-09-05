import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"period1_alignment/10077",
	"period1_physics_margin/274199",
	"period1_physics/274244",
};

string dataset_labels[] = {
	"alignment run",
	"$15\un{\si} + 0.5\un{mm}$",
	"$15\un{\si}$",
};

string rp = "L_1_F";

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(dataset_labels[dsi], "_", "\_") + "}");
}

NewRow();

for (int dsi : datasets.keys)
{
	NewPad("$x\ung{mm}$", "$y\ung{mm}$");

	scale(Linear, Linear, Log);

	RootGetObject(topDir + datasets[dsi]+"/distributions.root", "after selection/h2_y_vs_x_" + rp + "_sel");
	draw(robj, "def");

	limits((0, -15), (20, +15), Crop);
}

GShipout(hSkip=1mm, vSkip=1mm);
