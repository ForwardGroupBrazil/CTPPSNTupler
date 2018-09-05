import root;
import pad_layout;
include "../fills_samples.asy";

string topDir = "../../";

//----------------------------------------------------------------------------------------------------

InitDataSets("period1");

string ref_label[];

ref_label.push("10077");
ref_label.push("10079");
ref_label.push("10081");

string methods[];
pen method_pens[];

//methods.push("method y"); method_pens.push(blue);
methods.push("method x"); method_pens.push(red);

int rp_ids[];
string rps[], rp_labels[];
real rp_shift_m[], rp_shift_no_m[];
rp_ids.push(3); rps.push("L_1_F"); rp_labels.push("L-210-fr-hr"); rp_shift_m.push(-3.65); rp_shift_no_m.push(-4.17);
rp_ids.push(2); rps.push("L_1_N"); rp_labels.push("L-210-nr-hr"); rp_shift_m.push(-1.15); rp_shift_no_m.push(-0.00);
rp_ids.push(102); rps.push("R_1_N"); rp_labels.push("R-210-nr-hr"); rp_shift_m.push(-3.32); rp_shift_no_m.push(-3.93);
rp_ids.push(103); rps.push("R_1_F"); rp_labels.push("R-210-fr-hr"); rp_shift_m.push(-2.96); rp_shift_no_m.push(-3.57);

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
		int rp_id = rp_ids[rpi];

		for (int dsi : fill_data[fdi].datasets.keys)
		{
			string dataset = fill_data[fdi].datasets[dsi].tag;
			int dataset_idx = fill_data[fdi].datasets[dsi].idx;

			write("        " + dataset);
	
			mark m = (find(dataset, "margin") != -1) ? mSq+4pt+false : mCi+3pt;
	
			for (int mi : methods.keys)
			{
				real S1=0, Ss=0, Su=0;
	
				for (int ri : ref_label.keys)
				{
					RootObject obj = RootGetObject(topDir + dataset+"/match.root",
						rps[rpi] + "/" + ref_label[ri] + "/" + methods[mi] + "/g_results", error = false);
	
					if (!obj.valid)
						continue;
	
					real ax[] = { 0. };
					real ay[] = { 0. };
					obj.vExec("GetPoint", 0, ax, ay); real bsh = ay[0];
					obj.vExec("GetPoint", 1, ax, ay); real bsh_unc = ay[0];
		
					S1 += 1;
					Ss += bsh;		
					Su += bsh_unc;		
				}
				
				real x = fdi;

				real m_sh = Ss / S1;
				real u_sh = Su / S1;
	
				bool pointValid = (fabs(m_sh) > 0.01);
	
				pen p = StdPen(dataset_idx+1);
	
				if (pointValid)
				{
					draw((x, m_sh), m + p);
					draw((x, m_sh-u_sh)--(x, m_sh+u_sh), p);
				}
			}
		}
	}

	xlimits(-1, fill_data.length, Crop);

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
