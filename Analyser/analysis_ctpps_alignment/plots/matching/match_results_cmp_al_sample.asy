import root;
import pad_layout;
include "../fills_samples.asy";

string topDir = "../../";

//----------------------------------------------------------------------------------------------------

InitDataSets("period1_ps");

string ref_label[];
pen ref_pen[];
real ref_offset[];

ref_label.push("10077"); ref_pen.push(blue); ref_offset.push(-0.2);
ref_label.push("10079"); ref_pen.push(red); ref_offset.push(-0.0);
ref_label.push("10081"); ref_pen.push(heavygreen); ref_offset.push(+0.2);

string methods[];
mark method_markers[];

//methods.push("method y"); method_markers.push(mSq+2pt);
methods.push("method x"); method_markers.push(mCi+3pt);

string rps[], rp_labels[];
real rp_shift_m[], rp_shift_no_m[];
rps.push("L_1_F"); rp_labels.push("L-210-fr-hr"); rp_shift_m.push(-3.65); rp_shift_no_m.push(-4.17);
rps.push("L_1_N"); rp_labels.push("L-210-nr-hr"); rp_shift_m.push(-1.15); rp_shift_no_m.push(-0.00);
rps.push("R_1_N"); rp_labels.push("R-210-nr-hr"); rp_shift_m.push(-3.32); rp_shift_no_m.push(-3.93);
rps.push("R_1_F"); rp_labels.push("R-210-fr-hr"); rp_shift_m.push(-2.96); rp_shift_no_m.push(-3.57);

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

for (int rpi : rps.keys)
{
	write(rps[rpi]);

	if (rpi == 2)
		NewRow();

	NewPad("fill", "horizontal shift$\ung{mm}$");

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

	for (int fdi : fill_data.keys)
	{
		write(format("    %i", fill_data[fdi].fill));

		int fill = fill_data[fdi].fill; 

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			write("        " + dataset);

			mark m = (find(dataset, "margin") != -1) ? mSq+4pt+false : mCi+3pt;
	
			for (int ri : ref_label.keys)
			{
				for (int mi : methods.keys)
				{
					RootObject obj = RootGetObject(topDir + dataset+"/match.root",
						rps[rpi] + "/" + ref_label[ri] + "/" + methods[mi] + "/g_results", error=false);
					if (!obj.valid)
						continue;
	
					real ax[] = { 0. };
					real ay[] = { 0. };
					obj.vExec("GetPoint", 0, ax, ay); real bsh = ay[0];
					obj.vExec("GetPoint", 1, ax, ay); real bsh_unc = ay[0];
		
					real x = fdi + ref_offset[ri];
		
					pen p = ref_pen[ri];
	
					bool pointValid = (fabs(bsh) > 1);
				
					if (pointValid)
					{
						draw((x, bsh), m + p);
						draw((x, bsh-bsh_unc)--(x, bsh+bsh_unc), p);
					}
				}
			}
		}
	}

	xlimits(-1, fill_data.length);

	AttachLegend("{\SetFontSizesXX " + rp_labels[rpi] + "}");
}

//----------------------------------------------------------------------------------------------------

NewPad(false);

AddToLegend("run with margin", mSq+4pt+false);
AddToLegend("run without margin", mCi+3pt);

AddToLegend("<{\it TOTEM run in}");
AddToLegend("<{\it calibration fill}");
for (int ri : ref_label.keys)
	AddToLegend(ref_label[ri], ref_pen[ri]);

AttachLegend();

GShipout(hSkip=5mm, vSkip=1mm);
