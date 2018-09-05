import root;
import pad_layout;

string topDir = "../../";

string dataset = "period1_alignment/10077";

string selections[] = {
	"before selection",
	"after selection"
};

string selections_suffix[] = {
	"_no_sel",
	"_sel"
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

for (int seli : selections.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX " + replace(selections[seli], "_", "\_") + "}");

	for (int rpi : rps.keys)
	{
		NewPad("$x\ung{mm}$", "$y\ung{mm}$");

		scale(Linear, Linear, Log);

		RootGetObject(topDir + dataset+"/distributions.root", selections[seli] + "/h2_y_vs_x_" + rps[rpi] + selections_suffix[seli]);
		draw(robj, "def");

		limits((0, -15), (20, +15), Crop);
	}
}

GShipout(hSkip=1mm, vSkip=1mm);
