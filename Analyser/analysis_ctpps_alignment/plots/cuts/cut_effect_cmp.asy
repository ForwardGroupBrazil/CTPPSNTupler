import root;
import pad_layout;

string topDir = "../../";

string dataset = "period1_physics/fill_4985";

string rp = "L_1_N";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

string f = topDir + "/" + dataset + "/distributions.root";

NewPad("$x\ung{mm}$", "$y\ung{mm}$");
scale(Linear, Linear, Log);
draw(RootGetObject(f, "before selection/h2_y_vs_x_"+rp+"_no_sel"));
limits((5, -15), (20, +15), Crop);
AttachLegend("before");

NewPad("$x\ung{mm}$", "$y\ung{mm}$");
scale(Linear, Linear, Log);
draw(RootGetObject(f, "after selection/h2_y_vs_x_"+rp+"_sel"));
limits((5, -15), (20, +15), Crop);
AttachLegend("after");
