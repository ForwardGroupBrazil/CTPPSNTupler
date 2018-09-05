import root;
import pad_layout;
include "../fills_samples.asy";

string topDir = "../../";

//----------------------------------------------------------------------------------------------------

InitDataSets("period1_ps");
//InitDataSets("period2_ps");

string reco_dir = "none/without cuts";

string projections[];
pen pr_pens[];
projections.push("x"); pr_pens.push(blue);
projections.push("y"); pr_pens.push(red);

int rp_ids[];
string rps[], rp_labels[];
rp_ids.push(3); rps.push("L_1_F"); rp_labels.push("L-210-fr-hr");
rp_ids.push(2); rps.push("L_1_N"); rp_labels.push("L-210-nr-hr");
rp_ids.push(102); rps.push("R_1_N"); rp_labels.push("R-210-nr-hr");
rp_ids.push(103); rps.push("R_1_F"); rp_labels.push("R-210-fr-hr");

yTicksDef = RightTicks(1., 0.5);

xSizeDef = 15cm;

//----------------------------------------------------------------------------------------------------

string TickLabels(real x)
{
	if (x >=0 && x < fill_data.length)
	{
		return format("%i", fill_data[(int) x].fill);
	} else {
		return "";
	}
}

xTicksDef = LeftTicks(rotate(90)*Label(""), TickLabels, Step=1, step=0);

//----------------------------------------------------------------------------------------------------

for (int rpi : rps.keys)
{
	write(rps[rpi]);

	if (rpi == 2)
		NewRow();

	NewPad("fill", "local track angle$\ung{mrad}$");

	for (int fdi : fill_data.keys)
	{
		write(format("    %i", fill_data[fdi].fill));

		int fill = fill_data[fdi].fill; 
		int rp_id = rp_ids[rpi];

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			write("        " + dataset);

			for (int pri : projections.keys)
			{
				RootObject obj = RootGetObject(topDir + dataset+"/reconstruction_test.root",
					reco_dir + format("/h_th_" + projections[pri] + "_%i", rp_ids[rpi]), error=false);
				if (!obj.valid)
					continue;

				real x = fdi;

				real mean = obj.rExec("GetMean") * 1e3;
				real rms = obj.rExec("GetRMS") * 1e3;
				real n = obj.rExec("GetEntries");

				if (n < 1)
					continue;

				real mean_unc = rms / sqrt(n);

				if (abs(mean) < 0.01)
					continue;

				pen p = pr_pens[pri];

				draw((x, mean), mCi+2pt+p);
				draw((x, mean-mean_unc)--(x, mean+mean_unc), p);
			}
		}
	}

	xlimits(-1, fill_data.length, Crop);

	AttachLegend("{\SetFontSizesXX " + rp_labels[rpi] + "}", E, E);
}

//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int pri : projections.keys)
{
	pen p = pr_pens[pri];
	AddToLegend(projections[pri], mCi+2pt+p);
}

AttachLegend();

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=5mm, vSkip=1mm);
