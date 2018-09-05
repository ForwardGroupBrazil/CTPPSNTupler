import root;
import pad_layout;
include "../fills_samples.asy";
include "../io_alignment_format.asy";

string topDir = "../../";

//----------------------------------------------------------------------------------------------------

InitDataSets("period1");

string methods[];
pen method_pens[];

methods.push("y alignment"); method_pens.push(red);

int rp_ids[];
string rps[], rp_labels[];
real rp_y_cen[];
rp_ids.push(3); rps.push("L_1_F"); rp_labels.push("L-210-fr-hr"); rp_y_cen.push(-0.5);
rp_ids.push(2); rps.push("L_1_N"); rp_labels.push("L-210-nr-hr"); rp_y_cen.push(+0.2);
rp_ids.push(102); rps.push("R_1_N"); rp_labels.push("R-210-nr-hr"); rp_y_cen.push(-0.3);
rp_ids.push(103); rps.push("R_1_F"); rp_labels.push("R-210-fr-hr"); rp_y_cen.push(-0.4);

yTicksDef = RightTicks(0.2, 0.1);

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

/*
for (int rpi : rps.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + replace(rps[rpi], "_", "\_") + "}");
}

NewRow();
*/

for (int rpi : rps.keys)
{
	write(rps[rpi]);

	if (rpi == 2)
		NewRow();

	NewPad("fill", "vertical shift$\ung{mm}$");

	if (rp_ids[rpi] != 2)
		currentpad.yTicks = RightTicks(0.05, 0.01);

	// TODO
	/*
	if (rp_shift_m[rpi] != 0)
	{
		real sh = rp_shift_m[rpi], unc = 0.15;
		real fill_min = -1, fill_max = 6;
		draw((fill_min, sh+unc)--(fill_max, sh+unc), black+dashed);
		draw((fill_min, sh)--(fill_max, sh), black+1pt);
		draw((fill_min, sh-unc)--(fill_max, sh-unc), black+dashed);
		draw((fill_max, sh-2*unc), invisible);
		draw((fill_max, sh+2*unc), invisible);
	}
	
	if (rp_shift_no_m[rpi] != 0)
	{
		real sh = rp_shift_no_m[rpi], unc = 0.15;
		real fill_min = 3, fill_max = fill_data.length;
		draw((fill_min, sh+unc)--(fill_max, sh+unc), black+dashed);
		draw((fill_min, sh)--(fill_max, sh), black+1pt);
		draw((fill_min, sh-unc)--(fill_max, sh-unc), black+dashed);
		//filldraw((fill_min, sh-unc)--(fill_max, sh-unc)--(fill_max, sh+unc)--(fill_min, sh+unc)--cycle, black+opacity(0.1), nullpen);
		draw((fill_max, sh-2*unc), invisible);
		draw((fill_max, sh+2*unc), invisible);
	}
	*/

	for (int fdi : fill_data.keys)
	{
		write(format("    %i", fill_data[fdi].fill));

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			write("        " + dataset);
	
			mark m = (find(dataset, "margin") != -1) ? mSq+4pt+false : mCi+3pt;

			AlignmentResults arc[];
			int ret = LoadAlignmentResults(topDir + dataset+"/y_alignment.out", arc);
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
	}

	if (rp_ids[rpi] == 2)
		limits((-1, -0.4), (fill_data.length, 1), Crop);
	else 
		limits((-1, rp_y_cen[rpi]-0.2), (fill_data.length, rp_y_cen[rpi]+0.2), Crop);

	//yaxis(XEquals(23.5, false), heavygreen);

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
