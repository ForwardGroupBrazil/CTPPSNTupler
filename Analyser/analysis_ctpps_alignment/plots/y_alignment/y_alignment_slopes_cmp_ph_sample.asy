import root;
import pad_layout;
include "../fills_samples.asy";

string topDir = "../../";

//----------------------------------------------------------------------------------------------------

InitDataSets("period1");

string ref_label[];

int rp_ids[];
string rps[], rp_labels[];
real rp_y_min[], rp_y_max[];
rp_ids.push(3); rps.push("L_1_F"); rp_labels.push("L-210-fr-hr"); rp_y_min.push(-0.04); rp_y_max.push(-0.01);
rp_ids.push(2); rps.push("L_1_N"); rp_labels.push("L-210-nr-hr"); rp_y_min.push(-0.07); rp_y_max.push(-0.02);
rp_ids.push(102); rps.push("R_1_N"); rp_labels.push("R-210-nr-hr"); rp_y_min.push(-0.03); rp_y_max.push(-0.00);
rp_ids.push(103); rps.push("R_1_F"); rp_labels.push("R-210-fr-hr"); rp_y_min.push(-0.03); rp_y_max.push(-0.00);

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

	NewPad("fill", "slope$\ung{rad}$");

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
	
			mark m = (find(dataset, "margin") != -1) ? mSq+4pt+false : mCi+3pt;
	
			RootObject obj = RootGetObject(topDir + dataset+"/y_alignment.root",
				rps[rpi] + "/p_y_vs_x|ff_pol1", error = false);
	
			if (!obj.valid)
				continue;
	
			real x = fdi;

			real y = obj.rExec("GetParameter", 1);
			real y_unc = obj.rExec("GetParError", 1);

			pen p = StdPen(dataset_idx + 1);

			{
				draw((x, y), m + p);
				draw((x, y-y_unc)--(x, y+y_unc), p);
			}
		}
	}

	limits((-1, rp_y_min[rpi]), (fill_data.length, rp_y_max[rpi]), Crop);

	AttachLegend("{\SetFontSizesXX " + rp_labels[rpi] + "}");
}

//----------------------------------------------------------------------------------------------------

NewPad(false);

AddToLegend("run with margin", mSq+4pt+false);
AddToLegend("run without margin", mCi+3pt);

AddToLegend("main sample", red);
AddToLegend("validation sample", blue);

AttachLegend();

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=5mm, vSkip=1mm);
