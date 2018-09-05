void ApplySettings()
{
	lhc_fill = 5048;

	cut1_a = -1.018; cut1_c = -0.751; cut1_si = 0.20;
	cut2_a = -0.986; cut2_c = +0.657; cut2_si = 0.20;
	cut3_a = -1.153; cut3_c = +0.352; cut3_si = 0.15;
	cut4_a = -1.139; cut4_c = +0.044; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(13.3, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(12.0, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.8, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(10.0, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(8.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(8.5, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.5, 10.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(5.0, 10.0);
}
