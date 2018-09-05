import root;
import pad_layout;

string topDir = "../../";

string dataset = "period1_physics/fill_4985";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

string f = topDir + "/" + dataset + "/distributions.root";

NewPad("$x(\hbox{L-210-nr-hr})\ung{mm}$", "$x(\hbox{L-210-fr-hr})\ung{mm}$");
scale(Linear, Linear, Log);
draw(RootGetObject(f, "cuts/cut 1/canvas_before|h2_x_L_1_F_vs_x_L_1_N"));
draw(RootGetObject(f, "cuts/cut 1/canvas_before|g_up"), "l", magenta+1pt);
draw(RootGetObject(f, "cuts/cut 1/canvas_before|g_down"), "l", magenta+1pt);
limits((5, 5), (20, 20), Crop);

NewPad("$x(\hbox{R-210-nr-hr})\ung{mm}$", "$x(\hbox{R-210-fr-hr})\ung{mm}$");
scale(Linear, Linear, Log);
draw(RootGetObject(f, "cuts/cut 2/canvas_before|h2_x_R_1_F_vs_x_R_1_N"));
draw(RootGetObject(f, "cuts/cut 2/canvas_before|g_up"), "l", magenta+1pt);
draw(RootGetObject(f, "cuts/cut 2/canvas_before|g_down"), "l", magenta+1pt);
limits((5, 5), (20, 20), Crop);
