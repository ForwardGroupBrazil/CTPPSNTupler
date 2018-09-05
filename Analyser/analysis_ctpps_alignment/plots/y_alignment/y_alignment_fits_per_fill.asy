import root;
import pad_layout;

string topDir = "../../";

include "../fills_samples.asy";

//InitDataSets("period1_ps");
InitDataSets("period1_ss");

int rp_ids[];
string rps[], rp_labels[];
real rp_x_min[], rp_x_max[], rp_y_cen[];
rp_ids.push(3); rps.push("L_1_F"); rp_labels.push("L-210-fr-hr"); rp_x_min.push(-3.65); rp_x_max.push(-4.17); rp_y_cen.push(-0.75);
rp_ids.push(2); rps.push("L_1_N"); rp_labels.push("L-210-nr-hr"); rp_x_min.push(-1.15); rp_x_max.push(-0.00); rp_y_cen.push(-0.4);
rp_ids.push(102); rps.push("R_1_N"); rp_labels.push("R-210-nr-hr"); rp_x_min.push(-3.32); rp_x_max.push(-3.93); rp_y_cen.push(-0.4);
rp_ids.push(103); rps.push("R_1_F"); rp_labels.push("R-210-fr-hr"); rp_x_min.push(-2.96); rp_x_max.push(-3.57); rp_y_cen.push(-0.6);

xSizeDef = 9cm;

xTicksDef = LeftTicks(1., 0.5);
yTicksDef = RightTicks(0.1, 0.05);

//----------------------------------------------------------------------------------------------------

for (int fdi : fill_data.keys)
{
	write(format("    %i", fill_data[fdi].fill));

	for (int dsi : fill_data[fdi].datasets.keys)
	{
		string dataset = fill_data[fdi].datasets[dsi].tag;
		int dataset_idx = fill_data[fdi].datasets[dsi].idx;
		
		NewRow();

		NewPad(false);
		label(replace(dataset, "_", "\_"));

		for (int rpi : rps.keys)
		{
			NewPad("$x\ung{mm}$", "mean of $y\ung{mm}$");
		
			RootObject profile = RootGetObject(topDir + dataset + "/y_alignment.root", rps[rpi] + "/p_y_vs_x");
			RootObject fit = RootGetObject(topDir + dataset + "/y_alignment.root", rps[rpi] + "/p_y_vs_x|ff_pol1", error=false);

			if (!fit.valid)
				continue;

			draw(profile, "eb", red);

			TF1_x_min = -inf;
			TF1_x_max = +inf;
			draw(fit, "def", blue+2pt);

			TF1_x_min = 0;
			TF1_x_max = 13;
			draw(fit, "def", blue+dashed);

			//real y_cen = rp_y_cen[rpi];
			real y_cen = fit.rExec("Eval", 8.);
		
			limits((0, y_cen - 0.3), (13, y_cen + 0.3), Crop);
		
			AttachLegend(rp_labels[rpi]);
		}
	}
}
