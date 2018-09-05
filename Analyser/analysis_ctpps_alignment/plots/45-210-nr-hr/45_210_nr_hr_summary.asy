import root;
import pad_layout;
include "../fills_samples.asy";
include "../io_alignment_format.asy";

string topDir = "../../";

//----------------------------------------------------------------------------------------------------

InitDataSets("period1");
//InitDataSets("period2_ps");

string methods[];
pen method_pens[];
//methods.push("method y"); method_pens.push(blue);
methods.push("method x"); method_pens.push(red);

string reco_dir = "none/without cuts";

string projections[];
pen pr_pens[];
projections.push("x"); pr_pens.push(heavygreen);
projections.push("y"); pr_pens.push(magenta);

int rp_ids[];
string rps[], rp_labels[];
//rp_ids.push(3); rps.push("L_1_F"); rp_labels.push("L-210-fr-hr");
rp_ids.push(2); rps.push("L_1_N"); rp_labels.push("L-210-nr-hr");
//rp_ids.push(102); rps.push("R_1_N"); rp_labels.push("R-210-nr-hr");
//rp_ids.push(103); rps.push("R_1_F"); rp_labels.push("R-210-fr-hr");

yTicksDef = RightTicks(1., 0.5);

xSizeDef = 25cm;
ySizeDef = 5cm;

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
//----------------------------------------------------------------------------------------------------

for (int rpi : rps.keys)
{
	//write(rps[rpi]);

	if (rpi == 2)
		NewRow();

	NewPad("fill", "horizontal shift$\ung{mm}$");
	
	currentpad.yTicks = RightTicks(0.5, 0.1);

	for (int fdi : fill_data.keys)
	{
		//write(format("    %i", fill_data[fdi].fill));

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			//write("        " + dataset);
	
			mark m = (find(dataset, "margin") != -1) ? mSq+4pt+false : mCi+3pt;

			AlignmentResults arc[];
			int ret = LoadAlignmentResults(topDir + dataset+"/process_alignments.out", arc);
			if (ret != 0)
				continue;
	
			for (int mi : methods.keys)
			{
				bool method_found = false;
				AlignmentResults ar;
				for (AlignmentResults ari : arc)
				{
					if (ari.label == methods[mi])
					{
						method_found = true;
						ar = ari;
						break;
					}
				}

				if (!method_found)
					continue;

				if (!ar.results.initialized(rp_ids[rpi]))
					continue;

				real v = ar.results[rp_ids[rpi]].sh_x;

				bool point_valid = (fabs(v) > 1);

				if (point_valid)
				{
					pen p = StdPen(dataset_idx + 1);
					draw((fdi, v), m + p);
				}
			}
		}
	}

	xlimits(-1, fill_data.length, Crop);

	AddToLegend("run with margin", mSq+4pt+false);
	AddToLegend("run without margin", mCi+3pt);

	AddToLegend("main sample", red);
	AddToLegend("validation sample", blue);

	AttachLegend("{" + rp_labels[rpi] + "}");
}

//----------------------------------------------------------------------------------------------------

NewRow();

for (int rpi : rps.keys)
{
	//write(rps[rpi]);

	if (rpi == 2)
		NewRow();

	NewPad("fill", "vertical shift$\ung{mm}$");
	
	if (rp_ids[rpi] != 2)
		currentpad.yTicks = RightTicks(0.05, 0.01);
	else
		currentpad.yTicks = RightTicks(0.2, 0.1);

	for (int fdi : fill_data.keys)
	{
		//write(format("    %i", fill_data[fdi].fill));

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			//write("        " + dataset);
	
			mark m = (find(dataset, "margin") != -1) ? mSq+4pt+false : mCi+3pt;

			AlignmentResults arc[];
			int ret = LoadAlignmentResults(topDir + dataset+"/y_alignment.out", arc);
			if (ret != 0)
				continue;
	
			bool method_found = false;
			AlignmentResults ar;
			for (AlignmentResults ari : arc)
			{
				if (ari.label == "y alignment")
				{
					method_found = true;
					ar = ari;
					break;
				}
			}

			if (!method_found)
				continue;

			if (!ar.results.initialized(rp_ids[rpi]))
				continue;

			real v = ar.results[rp_ids[rpi]].sh_y;
			real v_unc = ar.results[rp_ids[rpi]].sh_y_unc;

			bool point_valid = (fabs(v) != 0.);

			if (point_valid)
			{
				pen p = StdPen(dataset_idx + 1);

				draw((fdi, v), m + p);
				draw((fdi, v-v_unc)--(fdi, v+v_unc), p);
			}
		}
	}

	/*
	if (rp_ids[rpi] == 2)
		limits((-1, -0.4), (fill_data.length, 1), Crop);
	else 
		limits((-1, rp_y_cen[rpi]-0.2), (fill_data.length, rp_y_cen[rpi]+0.2), Crop);
	*/

	//yaxis(XEquals(23.5, false), heavygreen);
	
	xlimits(-1, fill_data.length, Crop);

	AttachLegend("{" + rp_labels[rpi] + "}");
}

//----------------------------------------------------------------------------------------------------

NewRow();

for (int rpi : rps.keys)
{
	//write(rps[rpi]);

	if (rpi == 2)
		NewRow();

	NewPad("fill", "slope of $\langle y\rangle$ vs.~$x$ $\ung{rad}$");

	currentpad.yTicks = RightTicks(0.01, 0.005);

	for (int fdi : fill_data.keys)
	{
		//write(format("    %i", fill_data[fdi].fill));

		int fill = fill_data[fdi].fill; 
		int rp_id = rp_ids[rpi];

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			//write("        " + dataset);
	
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

	//limits((-1, rp_y_min[rpi]), (fill_data.length, rp_y_max[rpi]), Crop);
	xlimits(-1, fill_data.length, Crop);

	AttachLegend("{" + rp_labels[rpi] + "}");
}

//----------------------------------------------------------------------------------------------------

NewRow();

for (int rpi : rps.keys)
{
	//write(rps[rpi]);

	if (rpi == 2)
		NewRow();

	NewPad("fill", "mean local track angle$\ung{mrad}$");

	for (int fdi : fill_data.keys)
	{
		//write(format("    %i", fill_data[fdi].fill));

		int fill = fill_data[fdi].fill; 
		int rp_id = rp_ids[rpi];

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			//write("        " + dataset);

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

				draw((x, mean), mCi+3pt+p);
				draw((x, mean-mean_unc)--(x, mean+mean_unc), p);
			}
		}
	}

	xlimits(-1, fill_data.length, Crop);

	for (int pri : projections.keys)
	{
		AddToLegend(projections[pri], mCi+3pt+pr_pens[pri]);
	}

	AttachLegend("{" + rp_labels[rpi] + "}", E, E);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

GShipout(hSkip=5mm, vSkip=1mm);
